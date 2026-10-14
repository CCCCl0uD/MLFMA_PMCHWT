// RunLogger.h
#pragma once
#ifndef RUNLOGGER_H
#define RUNLOGGER_H

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <streambuf>
#include <string>

class FilteredTeeBuf : public std::streambuf {
public:
	FilteredTeeBuf(std::streambuf* consoleBuf, std::ostream* logStream, std::mutex* logMutex)
		: consoleBuf_(consoleBuf), logStream_(logStream), logMutex_(logMutex) {
	}

	~FilteredTeeBuf() override {
		if (!currentLine_.empty()) {
			flushLineToLog();
		}
	}

protected:
	int overflow(int ch) override {
		if (ch == traits_type::eof()) {
			return traits_type::not_eof(ch);
		}

		const char c = static_cast<char>(ch);
		if (consoleBuf_) {
			consoleBuf_->sputc(c);
		}
		consumeForLog(c);
		return ch;
	}

	std::streamsize xsputn(const char* s, std::streamsize count) override {
		if (consoleBuf_) {
			consoleBuf_->sputn(s, count);
		}
		for (std::streamsize i = 0; i < count; ++i) {
			consumeForLog(s[i]);
		}
		return count;
	}

	int sync() override {
		if (consoleBuf_) {
			consoleBuf_->pubsync();
		}
		if (logStream_) {
			std::lock_guard<std::mutex> lock(*logMutex_);
			logStream_->flush();
		}
		return 0;
	}

private:
	std::streambuf* consoleBuf_ = nullptr;
	std::ostream* logStream_ = nullptr;
	std::mutex* logMutex_ = nullptr;
	std::string currentLine_;
	bool skipCurrentLine_ = false;

	static std::string trim(const std::string& s) {
		const auto begin = std::find_if_not(s.begin(), s.end(), [](unsigned char ch) {
			return std::isspace(ch) != 0;
			});
		const auto end = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char ch) {
			return std::isspace(ch) != 0;
			}).base();
		if (begin >= end) {
			return "";
		}
		return std::string(begin, end);
	}

	static bool contains(const std::string& text, const std::string& needle) {
		return text.find(needle) != std::string::npos;
	}

	static bool shouldFilterLine(const std::string& line) {
		const std::string t = trim(line);
		if (t.empty()) {
			return false;
		}

		const bool looksLikeProgressBar =
			t.front() == '[' &&
			contains(t, "]") &&
			contains(t, "%") &&
			(contains(t, "=>") || contains(t, ">"));
		if (looksLikeProgressBar) {
			return true;
		}

		if (contains(t, "Matrix solver:")) return true;
		if (contains(t, "MGMRES warning")) return true;
		if (contains(t, "GMRES converged")) return true;
		if (contains(t, "GMRES reached max iterations")) return true;
		if (contains(t, "GMRES did not converge")) return true;
		if (contains(t, "CGS converged")) return true;
		if (contains(t, "CGS reached max iterations")) return true;
		if (contains(t, "CGS breakdown")) return true;
		if (contains(t, "CGS did not converge")) return true;
		if (contains(t, "Final residual")) return true;
		if (contains(t, "final residual")) return true;
		if (contains(t, "exact residual")) return true;
		if (contains(t, "Restart cycle ")) return true;
		if (contains(t, "Cycle ") && contains(t, " Res:")) return true;
		if (contains(t, "Cycle\t") && contains(t, " Res:")) return true;
		if (contains(t, "Iter ") && contains(t, " Res:")) return true;
		if (contains(t, "Iter\t") && contains(t, " Res:")) return true;
		if (contains(t, "MLFMA-Moving time point") && contains(t, "RCS computed")) return true;
		return false;
	}

	void consumeForLog(char c) {
		if (!logStream_ || !logStream_->good()) {
			return;
		}

		if (c == '\r') {
			currentLine_.clear();
			skipCurrentLine_ = true;
			return;
		}
		if (c == '\n') {
			flushLineToLog();
			currentLine_.clear();
			skipCurrentLine_ = false;
			return;
		}
		currentLine_.push_back(c);
	}

	void flushLineToLog() {
		if (skipCurrentLine_ || shouldFilterLine(currentLine_)) {
			return;
		}
		std::lock_guard<std::mutex> lock(*logMutex_);
		(*logStream_) << currentLine_ << '\n';
	}
};

class ScopedRunLogger {
public:
	explicit ScopedRunLogger(const std::string& logFilePath)
		: logFilePath_(logFilePath) {
		namespace fs = std::filesystem;
		const fs::path path(logFilePath_);
		if (path.has_parent_path()) {
			fs::create_directories(path.parent_path());
		}

		logFile_.open(logFilePath_, std::ios::out);
		if (!logFile_.is_open()) {
			std::cerr << "Warning: cannot open run log file: " << logFilePath_ << std::endl;
			return;
		}

		oldCoutBuf_ = std::cout.rdbuf();
		oldCerrBuf_ = std::cerr.rdbuf();
		coutBuf_ = new FilteredTeeBuf(oldCoutBuf_, &logFile_, &logMutex_);
		cerrBuf_ = new FilteredTeeBuf(oldCerrBuf_, &logFile_, &logMutex_);
		std::cout.rdbuf(coutBuf_);
		std::cerr.rdbuf(cerrBuf_);
		enabled_ = true;
	}

	ScopedRunLogger(const ScopedRunLogger&) = delete;
	ScopedRunLogger& operator=(const ScopedRunLogger&) = delete;

	~ScopedRunLogger() {
		if (enabled_) {
			std::cout.flush();
			std::cerr.flush();
			std::cout.rdbuf(oldCoutBuf_);
			std::cerr.rdbuf(oldCerrBuf_);
		}
		delete coutBuf_;
		delete cerrBuf_;
	}

	const std::string& path() const {
		return logFilePath_;
	}

	bool enabled() const {
		return enabled_;
	}

private:
	std::string logFilePath_;
	std::ofstream logFile_;
	std::mutex logMutex_;
	std::streambuf* oldCoutBuf_ = nullptr;
	std::streambuf* oldCerrBuf_ = nullptr;
	FilteredTeeBuf* coutBuf_ = nullptr;
	FilteredTeeBuf* cerrBuf_ = nullptr;
	bool enabled_ = false;
};

#endif // RUNLOGGER_H
