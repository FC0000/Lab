module;

#include <cassert>

export module lab:sample;

import :core;
import :estimate;

/* internal */ namespace lab
{
	template<typename S>
	struct _sample_size_t { using type = size_t; };

	template<std::ranges::sized_range S>
	struct _sample_size_t<S> { using type = std::ranges::range_size_t<S>; };
}

export namespace lab
{
	template<stdr::range Sample>
	using sample_value_t = stdr::range_value_t<Sample>;

	template<stdr::range Sample>
	using sample_size_t = _sample_size_t<Sample>::type;

	namespace _detail
	{
		struct analyze_sample_cpo;

		template<typename ValueType, typename SizeType>
		struct sample2_result_t;

		template<typename ValueType, typename SizeType>
		struct sample1_result_t
		{
			using value_type = ValueType;
			using size_type = SizeType;

			friend struct analyze_sample_cpo;
			friend struct sample2_result_t<value_type, size_type>;
		private:
			size_type
				_size = 0;
			value_type
				_min = std::numeric_limits<value_type>::infinity(),
				_max = -_min,
				_mean = 0,
				_m2 = 0,
				_m3 = 0,
				_m4 = 0;

			void _update(value_type x)
			{
				_min = std::min(_min, x);
				_max = std::max(_max, x);

				++_size;
				value_type
					fp_size = static_cast<value_type>(_size),
					delta = x - _mean,
					r_delta = delta / fp_size,
					r_delta2 = r_delta * r_delta,
					term = delta * r_delta * (_size - 1);

				_mean += r_delta;
				_m4 += term * r_delta2 * (fp_size * fp_size - value_type(3 * _size + 3)) + 6 * r_delta2 * _m2 - 4 * r_delta * _m3;
				_m3 += term * r_delta * (_size - 2) - 3 * r_delta * _m2;
				_m2 += term;
			}
		public:
			value_type min() const { return _min; }
			value_type max() const { return _max; }
			size_type size() const
			{
				return _size;
			}
			value_type mean() const
			{
				assert(size() != 0);
				return _mean;
			}
			value_type variance() const
			{
				assert(size() > 1);
				return _m2 / (size() - 1);
			}
			value_type skewness() const
			{
				assert(size() > 2);
				value_type term = std::sqrt(value_type(size() - 1)) * _m3 / std::pow(_m2, 1.5);
				return term + 2 * term / (size() - 2);
			}
			// relative
			value_type kurtosis() const
			{
				assert(size() > 3);
				value_type term1 = _m4 / (_m2 * _m2) + 6 / value_type(size() + 1) - 3,
					fp_size = static_cast<value_type>(size()),
					fp_size_2 = fp_size * fp_size,
					term2 = fp_size_2 - value_type(5 * size() + 6);
				return fp_size_2 * term1 / term2 - term1 / term2 + 3;
			}
			value_type mean_variance() const
			{
				return variance() / size();
			}

			value_type mean_stderr() const
			{
				return std::sqrt(mean_variance());
			}

			estimate_t<value_type> mean_estimate() const
			{
				return estimate_t<value_type>(mean(), mean_variance());
			}

			value_type stddev() const
			{
				return std::sqrt(variance());
			}
		};

		template<typename ValueType, typename SizeType>
		struct sample2_result_t
		{
			using value_type = ValueType;
			using size_type = SizeType;

			friend struct analyze_sample_cpo;

			sample1_result_t<value_type, size_type> x, y;
		private:
			value_type _co_m = 0;

			void _update(std::pair<value_type, value_type> p)
			{
				value_type delta_x = p.first - x._mean;

				x._update(p.first);
				y._update(p.second);

				_co_m += delta_x * (p.second - y._mean);
			}
		public:		
			size_type size() const
			{
				return x.size();
			}
			value_type covariance() const
			{
				return _co_m / (size() - 1);
			}
		};
	}
}

/* internal */ namespace lab::_detail
{
	template<typename T>            inline constexpr bool _is_valid_pair                  = false;
	template<std::floating_point T> inline constexpr bool _is_valid_pair<std::pair<T, T>> = true;
}

export namespace lab
{
	namespace _detail
	{
		struct analyze_sample_cpo
		{
			template<typename Sample> requires std::floating_point<sample_value_t<Sample>>
			auto operator()(Sample&& sample) const
			{
				using value_type = sample_value_t<Sample>;
				using size_type = sample_size_t<Sample>;
				sample1_result_t<value_type, size_type> result;

				for (auto x : sample)
					result._update(x);

				return result;
			}

			template<typename Sample> requires _is_valid_pair<sample_value_t<Sample>>
			auto operator()(Sample&& sample) const
			{
				using value_type = sample_value_t<Sample>::first_type;
				using size_type = sample_size_t<Sample>;
				sample2_result_t<value_type, size_type> result;

				for (auto x : sample)
					result._update(x);

				return result;
			}
		};
	}
	
	inline constexpr _detail::analyze_sample_cpo analyze_sample;
}