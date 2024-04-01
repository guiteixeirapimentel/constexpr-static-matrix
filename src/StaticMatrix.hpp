#pragma once
#include <array>
#include <cstdint>

template <typename T, size_t nlines = 3, size_t ncols = 3>
class StaticMatrix
{
public:
    static constexpr size_t nlines_ = nlines;
    static constexpr size_t ncols_ = ncols;
    using InternalType = T;
    using ThisMatrixType = StaticMatrix<T, nlines, ncols>;

public:
    constexpr StaticMatrix() noexcept = default;
    constexpr ~StaticMatrix() noexcept = default;

    constexpr explicit StaticMatrix(const std::array<T, nlines * ncols> &initial) noexcept
        : data(initial)
    {
    }

    template <size_t i, size_t j>
    [[nodiscard]] inline constexpr T getVal() const noexcept
    {
        static_assert(i < nlines_);
        static_assert(j < ncols_);
        return data[j + (i * ncols_)];
    }

    template <size_t i, size_t j>
    constexpr inline void setVal(const T val) noexcept
    {
        static_assert(i < nlines_);
        static_assert(j < ncols_);
        data[j + (i * ncols_)] = val;
    }

    [[nodiscard]] inline constexpr T getVal(size_t i, size_t j) const noexcept
    {
        return data[j + (i * ncols_)];
    }

    constexpr inline void setVal(size_t i, size_t j, const T val) noexcept
    {
        data[j + (i * ncols_)] = val;
    }

    constexpr inline T& operator()(size_t i, size_t j) noexcept
    {
        return data[j + (i * ncols_)];
    }

    constexpr inline const T& operator()(size_t i, size_t j) const
    {
        return data[j + (i * ncols_)];
    }

    constexpr inline ThisMatrixType operator*(ThisMatrixType::InternalType scalar) noexcept
    {
        ThisMatrixType res;

        for(size_t i = 0; i < res.data.size(); i++)
        {
            res.data[i] = data[i] * scalar;
        }

        return res;
    }

    template <size_t rhs_nlines, size_t rhs_ncols>
    constexpr StaticMatrix<T, nlines, rhs_ncols> operator*(const StaticMatrix<T, rhs_nlines, rhs_ncols> &rhs) const noexcept
    {
        static_assert(ncols_ == rhs.nlines_);

        StaticMatrix<T, nlines_, rhs.ncols_> res;

        for (size_t i = 0; i < nlines_; i++)
        {
            for (size_t j = 0; j < rhs_ncols; j++)
            {
                T val = 0.0;

                for (size_t index = 0; index < rhs_nlines; index++)
                {
                    val += getVal(i, index) * rhs.getVal(index, j);
                }

                res.setVal(i, j, val);
            }
        }

        return res;
    }

    constexpr inline ThisMatrixType &operator+=(const ThisMatrixType &rhs) noexcept
    {
        for (size_t i = 0; i < data.size(); i++)
        {
            data[i] += rhs.data[i];
        }
        return *this;
    }

    constexpr inline ThisMatrixType operator+(const ThisMatrixType &rhs) const noexcept
    {
        ThisMatrixType res = *this;
        for (size_t i = 0; i < data.size(); i++)
        {
            res.data[i] += rhs.data[i];
        }
        return res;
    }

    constexpr inline ThisMatrixType &operator-=(const ThisMatrixType &rhs) noexcept
    {
        for (size_t i = 0; i < data.size(); i++)
        {
            data[i] -= rhs.data[i];
        }
        return *this;
    }

    constexpr inline ThisMatrixType operator-(const ThisMatrixType &rhs) const noexcept
    {
        ThisMatrixType res = *this;
        for (size_t i = 0; i < data.size(); i++)
        {
            res.data[i] -= rhs.data[i];
        }
        return res;
    }

    constexpr inline bool operator==(const ThisMatrixType &rhs) const noexcept
    {
        return data == rhs.data;
    }

    constexpr inline StaticMatrix<T, ncols_, nlines_> getTransposed() const noexcept
    {
        StaticMatrix<T, ncols_, nlines_> res{};

        for (size_t i = 0; i < nlines_; i++)
        {
            for (size_t j = 0; j < ncols_; j++)
            {
                res.setVal(j, i, getVal(i, j));
            }
        }

        return res;
    }

    static constexpr ThisMatrixType zeros()
    {
        ThisMatrixType res{};

        for (size_t i = 0; i < res.data.size(); i++)
        {
            res.data[i] = T{0};
        }

        return res;
    }

    static constexpr ThisMatrixType identity()
    {
        ThisMatrixType res{};

        for (size_t i = 0; i < nlines_; i++)
        {
            for (size_t j = 0; j < ncols_; j++)
            {
                res.setVal(i, j, i == j ? T{1} : T{0});
            }
        }

        return res;
    }

public:
    std::array<T, nlines_ * ncols_> data;
};