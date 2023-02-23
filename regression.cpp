module;

#include <cassert>

export module lab:regression;

import :core;
import :estimate;
import :sample;

export namespace lab
{
	namespace _detail
	{
		struct regression_cpo
		{
		private:
			template<typename ValueType, typename SizeType>
			struct regression_result_t
			{
				friend struct regression_cpo;

				using value_type = ValueType;
				using size_type = SizeType;
			private:
				value_type _slope, _intercept, _slope_stderr2, _intercept_stderr2, _R2;
				sample2_result_t<value_type, size_type> _sample;

				regression_result_t(value_type s, value_type i, value_type s_stderr2, value_type i_stderr2, value_type R2, sample2_result_t<value_type, size_type> const& sample)
					: _slope(s), _intercept(i), _slope_stderr2(s_stderr2), _intercept_stderr2(i_stderr2), _R2(R2), _sample(sample) {}
			public:
				value_type slope() const { return _slope; }
				value_type intercept() const { return _intercept; }
				value_type slope_stderr() const { return std::sqrt(_slope_stderr2); }
				value_type slope_stderr2() const { return _slope_stderr2; }
				estimate_t<value_type> slope_estimate() const { return { slope(), slope_stderr2() }; }
				value_type intercept_stderr() const { return std::sqrt(_intercept_stderr2); }
				value_type intercept_stderr2() const { return _intercept_stderr2; }
				estimate_t<value_type> intercept_estimate() const { return { intercept(), intercept_stderr2() }; }
				value_type R2() const { return _R2; }
				sample2_result_t<value_type, size_type> const& sample() const { return _sample;}
			};

			template<typename Sample>
			using _value_t = sample_value_t<Sample>::first_type;
		public:
			template<typename Sample>
			auto operator()(Sample const& sample, _value_t<Sample> y_variance = std::numeric_limits<_value_t<Sample>>::quiet_NaN()) const
			{
				using value_type = _value_t<Sample>;
				using size_type = sample_size_t<Sample>;

				auto s = analyze_sample(sample);

				auto size = s.size();
				assert(size > 2);

				value_type
					covariance = s.covariance(),
					x_mean = s.x.mean(),
					x_variance = s.x.variance(),
					y_mean = s.y.mean(),
					slope = covariance / x_variance,
					intercept = y_mean - slope * x_mean,
					y_post_variance = 0;

				for(auto [x, y] : sample)
				{
					value_type epsilon = y - (intercept + slope * x);
					y_post_variance += epsilon * epsilon;
				}
				y_post_variance /= size - 2;

				value_type
					y_var_or_post_var = std::isnan(y_variance) ? y_post_variance : y_variance,
					slope_stderr2 = y_var_or_post_var / ((size - 1) * x_variance),
					intercept_stderr2 = y_var_or_post_var / size + slope_stderr2 * x_mean * x_mean;

				return regression_result_t<value_type, size_type>(
					slope,
					intercept,
					slope_stderr2,
					intercept_stderr2,
					1 - y_post_variance / s.y.variance(),
					s);
			}
		};
	}

	inline constexpr _detail::regression_cpo regression;
}

