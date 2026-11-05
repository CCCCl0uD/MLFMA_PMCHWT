// NASReader.cpp
#include "NASReader.h"
#include "PathUtils.h"

#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

using namespace std;

namespace {
struct RawFaceRecord {
    int id = 0;
    int propertyId = 0;
    array<int, 3> vertices{};
    size_t sourceLine = 0;
};

string_view trimField(string_view value) {
    const size_t first = value.find_first_not_of(" \t\r");
    if (first == string_view::npos) return {};
    return value.substr(first, value.find_last_not_of(" \t\r") - first + 1);
}

string_view fixedField(string_view line, size_t offset, size_t width) {
    return offset < line.size() ? trimField(line.substr(offset, width)) : string_view{};
}

bool parseInteger(string_view field, int& value) {
    if (field.empty()) return false;
    if (field.front() == '+') field.remove_prefix(1);
    const auto parsed = from_chars(field.data(), field.data() + field.size(), value);
    return parsed.ec == errc{} && parsed.ptr == field.data() + field.size();
}

bool parseReal(string_view field, double& value) {
    if (field.empty()) return false;
    string number(field);
    for (char& c : number) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    if (number.find_first_of("Ee") == string::npos) {
        const size_t exponent = number.find_last_of("+-");
        if (exponent != string::npos && exponent != 0) number.insert(exponent, "E");
    }
    string_view normalized(number);
    if (normalized.front() == '+') normalized.remove_prefix(1);
    const auto parsed = from_chars(normalized.data(), normalized.data() + normalized.size(),
        value, chars_format::general);
    return parsed.ec == errc{} && parsed.ptr == normalized.data() + normalized.size()
        && isfinite(value);
}
} // namespace

Point CalculateTriangleCenter(const Point& a, const Point& b, const Point& c) {
    return {0, (a.x + b.x + c.x) / 3.0,
        (a.y + b.y + c.y) / 3.0, (a.z + b.z + c.z) / 3.0};
}

int readNasData(const string& filename, vector<Point>& points, vector<FaceElement>& faces) {
    MeshReadReport report;
    return readNasData(filename, points, faces, report);
}

int readNasData(const string& filename, vector<Point>& points,
    vector<FaceElement>& faces, MeshReadReport& report) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file: " << filename << '\n';
        return -1;
    }

    vector<Point> parsedPoints;
    vector<RawFaceRecord> rawFaces;
    unordered_map<int, size_t> nodeIndices;
    unordered_set<int> faceIds;
    string line;
    size_t lineNumber = 0;
    bool pendingGrid = false;
    Point pendingPoint{};
    string continuationLabel;

    const auto fail = [&](size_t sourceLine, const string& message) {
        cerr << filename << ':' << sourceLine << ": " << message << '\n';
        return -1;
    };

    while (getline(file, line)) {
        ++lineNumber;
        const size_t comment = line.find('$');
        if (comment != string::npos) line.resize(comment);
        const string_view text(line);
        if (trimField(text).empty()) continue;

        if (pendingGrid) {
            if (text.front() != '*') return fail(lineNumber, "Expected GRID* continuation");
            string_view label = fixedField(text, 0, 8);
            label = trimField(label.substr(1));
            if (!continuationLabel.empty() && label != continuationLabel)
                return fail(lineNumber, "GRID* continuation label mismatch");
            if (!parseReal(fixedField(text, 8, 16), pendingPoint.z))
                return fail(lineNumber, "Invalid GRID* Z coordinate");
            if (!nodeIndices.emplace(pendingPoint.id, parsedPoints.size()).second)
                return fail(lineNumber, "Duplicate GRID ID " + to_string(pendingPoint.id));
            parsedPoints.push_back(pendingPoint);
            pendingGrid = false;
            continue;
        }

        array<string_view, 10> fields{};
        const bool freeField = text.find(',') != string_view::npos;
        if (freeField) {
            size_t start = 0;
            for (auto& field : fields) {
                const size_t end = text.find(',', start);
                field = trimField(text.substr(start, end == string_view::npos ? end : end - start));
                if (end == string_view::npos) break;
                start = end + 1;
            }
        }
        else {
            fields[0] = fixedField(text, 0, 8);
            const size_t width = fields[0] == "GRID*" ? 16 : 8;
            const size_t count = fields[0] == "GRID*" ? 5 : fields.size();
            for (size_t i = 1; i < count; ++i) fields[i] = fixedField(text, 8 + (i - 1) * width, width);
        }

        if (fields[0] == "GRID*" || fields[0] == "GRID") {
            const bool largeField = fields[0] == "GRID*";
            if (largeField && freeField) return fail(lineNumber, "Free-field GRID* is not supported; use GRID");
            Point point{};
            int coordinateSystem = 0;
            if (!parseInteger(fields[1], point.id) || point.id <= 0 ||
                (!fields[2].empty() && !parseInteger(fields[2], coordinateSystem)) ||
                !parseReal(fields[3], point.x) || !parseReal(fields[4], point.y))
                return fail(lineNumber, "Invalid GRID ID, CP or coordinate");
            if (coordinateSystem != 0) return fail(lineNumber, "Nonzero GRID CP is not supported");
            if (largeField) {
                pendingPoint = point;
                pendingGrid = true;
                string_view label = fixedField(text, 72, 8);
                if (!label.empty() && label.front() == '*') label = trimField(label.substr(1));
                continuationLabel.assign(label);
            }
            else {
                if (!parseReal(fields[5], point.z)) return fail(lineNumber, "Invalid GRID Z coordinate");
                if (!nodeIndices.emplace(point.id, parsedPoints.size()).second)
                    return fail(lineNumber, "Duplicate GRID ID " + to_string(point.id));
                parsedPoints.push_back(point);
            }
        }
        else if (fields[0] == "CTRIA3") {
            RawFaceRecord face;
            face.sourceLine = lineNumber;
            if (!parseInteger(fields[1], face.id) || face.id <= 0 ||
                (!fields[2].empty() && (!parseInteger(fields[2], face.propertyId) || face.propertyId <= 0)))
                return fail(lineNumber, "Invalid CTRIA3 EID or PID");
            if (fields[2].empty()) face.propertyId = face.id;
            for (size_t j = 0; j < 3; ++j) {
                if (!parseInteger(fields[j + 3], face.vertices[j]) || face.vertices[j] <= 0)
                    return fail(lineNumber, "Invalid CTRIA3 node reference");
            }
            if (!faceIds.insert(face.id).second)
                return fail(lineNumber, "Duplicate CTRIA3 ID " + to_string(face.id));
            rawFaces.push_back(face);
        }
        else if (fields[0] == "INCLUDE" || fields[0] == "GRDSET" ||
            fields[0].substr(0, 4) == "CORD" || fields[0].substr(0, 5) == "CQUAD" ||
            fields[0].substr(0, 5) == "CTRIA" || fields[0] == "CTETRA" ||
            fields[0] == "CHEXA" || fields[0] == "CPENTA") {
            return fail(lineNumber, "Unsupported mesh card " + string(fields[0]));
        }
    }

    if (file.bad()) return fail(lineNumber, "NAS read error");
    if (pendingGrid) return fail(lineNumber, "Unfinished GRID* entry");
    if (parsedPoints.empty() || rawFaces.empty()) return fail(lineNumber, "No triangular mesh found");

    // No point is appended after this point. All face pointers borrow this storage.
    vector<FaceElement> parsedFaces;
    MeshReadReport parsedReport;
    parsedFaces.reserve(rawFaces.size());
    parsedReport.facePropertyIds.reserve(rawFaces.size());
    parsedReport.faceSourceLines.reserve(rawFaces.size());
    for (const auto& raw : rawFaces) {
        array<Point*, 3> vertices{};
        for (size_t j = 0; j < 3; ++j) {
            const auto found = nodeIndices.find(raw.vertices[j]);
            if (found == nodeIndices.end())
                return fail(raw.sourceLine, "Face " + to_string(raw.id) +
                    " references missing GRID " + to_string(raw.vertices[j]));
            vertices[j] = &parsedPoints[found->second];
        }
        if (vertices[0] == vertices[1] || vertices[1] == vertices[2] || vertices[2] == vertices[0])
            return fail(raw.sourceLine, "Repeated vertex in face " + to_string(raw.id));

        FaceElement face{};
        face.faceId = raw.id;
        face.vertex1 = vertices[0];
        face.vertex2 = vertices[1];
        face.vertex3 = vertices[2];
        const Point ab = *vertices[1] - *vertices[0];
        const Point ac = *vertices[2] - *vertices[0];
        Point normal{0, ab.y * ac.z - ab.z * ac.y,
            ab.z * ac.x - ab.x * ac.z, ab.x * ac.y - ab.y * ac.x};
        const double twiceArea = hypot(normal.x, normal.y, normal.z);
        face.triangularArea = twiceArea / 2.0;
        face.triangleCenter = CalculateTriangleCenter(*vertices[0], *vertices[1], *vertices[2]);
        if (!(face.triangularArea > 0.0) || !isfinite(twiceArea) ||
            !isfinite(face.triangleCenter.x) || !isfinite(face.triangleCenter.y) ||
            !isfinite(face.triangleCenter.z))
            return fail(raw.sourceLine, "Degenerate or non-finite face " + to_string(raw.id));
        normal = normal / twiceArea;
        if (abs(normal.x) < 1e-6) normal.x = 0;
        if (abs(normal.y) < 1e-6) normal.y = 0;
        if (abs(normal.z) < 1e-6) normal.z = 0;
        face.externalNormalVector = normal;
        parsedFaces.push_back(face);
        parsedReport.facePropertyIds.push_back(raw.propertyId);
        parsedReport.faceSourceLines.push_back(raw.sourceLine);
    }

    cout << "Info: NAS file read, points: " << parsedPoints.size()
        << ", faces: " << parsedFaces.size() << " (NASReader.cpp)\n";
    points.swap(parsedPoints);
    faces.swap(parsedFaces);
    report.swap(parsedReport);
    return 0;
}

bool writeOutputFiles(const string& nasPath, const vector<FaceElement>& faces, const vector<Point>& points) {
	string parentPath = getParentPath(nasPath);
	string baseName = getBaseName(nasPath);

	ofstream faceFile(parentPath + baseName + "_DEBUG_Faces.txt");
	if (!faceFile) {
		cerr << "Error creating faces output file" << endl;
		return false;
	}

	faceFile << "FaceID\tVertex1\tVertex2\tVertex3\tArea\tCenterX\tCenterY\tCenterZ\tExternalNormalVector(x)\tExternalNormalVector(y)\tExternalNormalVector(z)\n";
	for (const auto& face : faces) {
		faceFile << face.faceId << "\t"
			<< face.vertex1->id << "\t"
			<< face.vertex2->id << "\t"
			<< face.vertex3->id << "\t"
			<< face.triangularArea << "\t"
			<< face.triangleCenter.x << "\t"
			<< face.triangleCenter.y << "\t"
			<< face.triangleCenter.z << "\t"
			<< face.externalNormalVector.x << "\t"
			<< face.externalNormalVector.y << "\t"
			<< face.externalNormalVector.z << "\n";
	}

	cout << "Info: Face element data written to: " << parentPath + baseName + "_DEBUG_Faces.txt (NASReader.cpp)" << std::endl;

	ofstream pointFile(parentPath + baseName + "_DEBUG_Points.txt");
	if (!pointFile) {
		cerr << "Error creating points output file" << endl;
		return false;
	}

	pointFile << fixed << setprecision(6);
	pointFile << "PointID\tX\tY\tZ\n";
	for (const auto& point : points) {
		pointFile << point.id << "\t"
			<< point.x << "\t"
			<< point.y << "\t"
			<< point.z << "\n";
	}
	cout << "Info: Point data written to: " << parentPath + baseName + "_DEBUG_Points.txt (NASReader.cpp)" << std::endl;

	return true;
}