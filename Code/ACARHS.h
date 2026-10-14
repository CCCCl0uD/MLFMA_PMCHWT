// ACARHS.h
#pragma once
#ifndef ACA_RHS_H
#define ACA_RHS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ACARHS {
	struct Options {
		double tol = 1.0e-3;
		int maxRank = 0;
		double pivotEps = 1.0e-14;
	};

	struct Factors {
		int rows = 0;
		int cols = 0;
		int rank = 0;
		std::vector<std::complex<double>> U;
		std::vector<std::complex<double>> V;
	};

	inline double abs2(const std::complex<double>& z) {
		return std::norm(z);
	}

	inline bool containsIndex(const std::vector<int>& indices, int value) {
		return std::find(indices.begin(), indices.end(), value) != indices.end();
	}

	template<typename ComputeColumnFunc, typename ComputeEntryFunc>
	Factors factorize(
		int rows,
		int cols,
		const Options& opt,
		ComputeColumnFunc computeColumn,
		ComputeEntryFunc computeEntry)
	{
		if (rows <= 0 || cols <= 0) {
			throw std::invalid_argument("ACARHS::factorize requires positive matrix dimensions.");
		}
		if (opt.tol <= 0.0) {
			throw std::invalid_argument("ACARHS::factorize requires a positive tolerance.");
		}
		if (opt.pivotEps <= 0.0) {
			throw std::invalid_argument("ACARHS::factorize requires a positive pivot threshold.");
		}

		const int rankLimit = (opt.maxRank > 0) ?
			std::min({ opt.maxRank, rows, cols }) :
			std::min(rows, cols);

		Factors factors;
		factors.rows = rows;
		factors.cols = cols;
		factors.U.reserve(static_cast<size_t>(rankLimit) * rows);
		factors.V.reserve(static_cast<size_t>(rankLimit) * cols);

		int pivotRow = 0;
		double approxNorm = 0.0;
		std::vector<int> selectedRows;
		std::vector<int> selectedCols;
		selectedRows.push_back(pivotRow);

		for (int rank = 0; rank < rankLimit; ++rank) {
			std::vector<std::complex<double>> rowResidual(
				cols, std::complex<double>(0.0, 0.0));

			int pivotCol = -1;
			double pivotColAbs2 = 0.0;
			for (int col = 0; col < cols; ++col) {
				std::complex<double> value = computeEntry(pivotRow, col);
				for (int p = 0; p < factors.rank; ++p) {
					value -= factors.U[static_cast<size_t>(p) * rows + pivotRow] *
						factors.V[static_cast<size_t>(p) * cols + col];
				}
				rowResidual[col] = value;

				const double valueAbs2 = abs2(value);
				if (!containsIndex(selectedCols, col) && valueAbs2 > pivotColAbs2) {
					pivotColAbs2 = valueAbs2;
					pivotCol = col;
				}
			}

			if (pivotCol < 0 || pivotColAbs2 <= opt.pivotEps * opt.pivotEps) {
				std::cout << "ACA-RHS stopped: pivot column residual is below threshold at rank "
					<< factors.rank << std::endl;
				break;
			}
			selectedCols.push_back(pivotCol);

			const std::complex<double> pivot = rowResidual[pivotCol];
			std::vector<std::complex<double>> newV(cols);
			double norm2V = 0.0;
			for (int col = 0; col < cols; ++col) {
				newV[col] = rowResidual[col] / pivot;
				norm2V += abs2(newV[col]);
			}

			std::vector<std::complex<double>> column;
			computeColumn(pivotCol, column);
			if (static_cast<int>(column.size()) != rows) {
				throw std::runtime_error("ACARHS::factorize computeColumn returned an unexpected length.");
			}

			std::vector<std::complex<double>> newU(rows);
			double norm2U = 0.0;
			int nextPivotRow = -1;
			double nextPivotRowAbs2 = 0.0;
			for (int row = 0; row < rows; ++row) {
				std::complex<double> value = column[row];
				for (int p = 0; p < factors.rank; ++p) {
					value -= factors.U[static_cast<size_t>(p) * rows + row] *
						factors.V[static_cast<size_t>(p) * cols + pivotCol];
				}
				newU[row] = value;
				norm2U += abs2(value);

				const double valueAbs2 = abs2(value);
				if (!containsIndex(selectedRows, row) && valueAbs2 > nextPivotRowAbs2) {
					nextPivotRowAbs2 = valueAbs2;
					nextPivotRow = row;
				}
			}

			factors.U.insert(factors.U.end(), newU.begin(), newU.end());
			factors.V.insert(factors.V.end(), newV.begin(), newV.end());
			++factors.rank;

			const double newContribution = norm2U * norm2V;
			approxNorm += newContribution;
			std::cout << "ACA-RHS rank " << factors.rank
				<< ", pivot row = " << pivotRow
				<< ", pivot col = " << pivotCol
				<< ", contribution = " << newContribution << std::endl;

			if (factors.rank > 1 &&
				newContribution <= opt.tol * opt.tol * std::max(approxNorm, std::numeric_limits<double>::min())) {
				std::cout << "ACA-RHS converged at rank " << factors.rank << std::endl;
				break;
			}
			if (nextPivotRow < 0 || nextPivotRowAbs2 <= opt.pivotEps * opt.pivotEps) {
				std::cout << "ACA-RHS stopped: pivot row residual is below threshold at rank "
					<< factors.rank << std::endl;
				break;
			}

			pivotRow = nextPivotRow;
			selectedRows.push_back(pivotRow);
		}

		if (factors.rank <= 0) {
			throw std::runtime_error("ACA-RHS produced zero rank factors.");
		}
		return factors;
	}

	inline std::vector<std::complex<double>> reconstructColumn(
		const Factors& factors,
		int col)
	{
		if (col < 0 || col >= factors.cols) {
			throw std::out_of_range("ACARHS::reconstructColumn column index is out of range.");
		}

		std::vector<std::complex<double>> column(
			factors.rows, std::complex<double>(0.0, 0.0));
		for (int m = 0; m < factors.rank; ++m) {
			const std::complex<double> coeff = factors.V[static_cast<size_t>(m) * factors.cols + col];
			for (int row = 0; row < factors.rows; ++row) {
				column[row] += factors.U[static_cast<size_t>(m) * factors.rows + row] * coeff;
			}
		}
		return column;
	}
}

#endif // ACA_RHS_H
