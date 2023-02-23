module;

#include <cassert>

export module lab:estimate;

import :core;

export namespace lab
{
	template<typename T = double>
	struct estimate_t
	{
		using value_type = T;
		static_assert(std::floating_point<value_type>);
	private:
		value_type _value, _variance;
	public:
		constexpr estimate_t()
			: _value(std::numeric_limits<value_type>::quiet_NaN()), _variance(_value)
		{}

		constexpr estimate_t(value_type value, value_type variance)
			: _value(value), _variance(variance)
		{
			assert(variance >= 0);
		}

		constexpr value_type value() const { return _value; }
		constexpr value_type variance() const { return _variance; }
		constexpr value_type stddev() const { return std::sqrt(variance()); }
	};
}

export namespace std
{
	template<typename T, typename CharT>
	struct formatter<lab::estimate_t<T>, CharT> : formatter<T, CharT>
	{
		template<typename FormatContext>
		constexpr auto format(lab::estimate_t<T> estimate, FormatContext& format_context) const
		{
			auto out = formatter<T, CharT>::format(estimate.value(), format_context);
			out = ranges::copy(" +- ", out).out;
			format_context.advance_to(out);
			return formatter<T, CharT>::format(estimate.stddev(), format_context);
		}
	};
}

export namespace lab
{
	template<typename value_type, size_t N, typename Range = decltype(stdv::repeat(0))>
	auto estimate(auto const& function, std::array<value_type, N> const& arguments, std::array<value_type, N> const& true_values, Range&& covariance_triangular_matrix = stdv::repeat(0))
	{
		static_assert(stdr::range<Range>);
	
		auto [value, derivative] = [&]<size_t... In>(std::index_sequence<In...>)
		{
			return std::pair(function.value_at(arguments[In]...), function.derivative_at(true_values[In]...));
		}(std::make_index_sequence<N>());

		value_type variance = 0, compensation = 0;
		auto compensated_add = [&](value_type new_term)
		{
			value_type t = variance + new_term;
			if (std::abs(variance) >= std::abs(new_term))
				compensation += (variance - t) + new_term;
			else
				compensation += (new_term - t) + variance;
			variance = t;
		};
		auto cov_it = stdr::begin(covariance_triangular_matrix);
		for (size_t i = 0; i != N; ++i)
		{
			compensated_add(derivative[i] * derivative[i] * (*cov_it++));
			for (size_t j = i; j != N; ++j)
				compensated_add(derivative[i] * derivative[j] * (*cov_it++));
		}
		variance += compensation;
		return estimate_t(value, variance);
	}

	template<typename value_type, size_t N, typename Range = decltype(stdv::repeat(0))>
	auto estimate(auto const& function, std::array<value_type, N> const& arguments, Range&& covariance_triangular_matrix = stdv::repeat(0))
	{
		return estimate(function, arguments, arguments, covariance_triangular_matrix);
	}


	template<typename value_type, size_t N, typename Range = decltype(stdv::repeat(0))>
	auto estimate(auto const& function, std::array<estimate_t<value_type>, N> const& arguments, std::array<value_type, N> const& true_values, Range&& covariance_strictly_triangular_matrix = stdv::repeat(0))
	{
		static_assert(stdr::range<Range>);

		auto [value, derivative] = [&]<size_t... In>(std::index_sequence<In...>)
		{
			return std::pair(function.value_at(arguments[In].value()...), function.derivative_at(true_values[In]...));
		}(std::make_index_sequence<N>());

		value_type variance = 0, compensation = 0;
		auto compensated_add = [&](value_type new_term)
		{
			value_type t = variance + new_term;
			if (std::abs(variance) >= std::abs(new_term))
				compensation += (variance - t) + new_term;
			else
				compensation += (new_term - t) + variance;
			variance = t;
		};
		auto cov_it = stdr::begin(covariance_strictly_triangular_matrix);
		for (size_t i = 0; i != N; ++i)
		{
			compensated_add(derivative[i] * derivative[i] * arguments[i].variance());
			for (size_t j = i; j != N; ++j)
				compensated_add(2 * derivative[i] * derivative[j] * (*cov_it++));
		}
		variance += compensation;
		return estimate_t(value, variance);
	}
	template<typename value_type, size_t N, typename Range = decltype(stdv::repeat(0))>
	auto estimate(auto const& function, std::array<estimate_t<value_type>, N> const& arguments, Range&& covariance_strictly_triangular_matrix = stdv::repeat(0))
	{
		std::array<value_type, N> true_values;
		stdr::copy(arguments | stdv::transform([](auto x) {return x.value(); }), true_values.begin());
		return estimate(function, arguments, true_values, covariance_strictly_triangular_matrix);
	}
}