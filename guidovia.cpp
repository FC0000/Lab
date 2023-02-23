/*
	Esperienza guidovia
*/

import lab;

namespace stdv = std::views;
namespace stdr = std::ranges;
namespace stdf = std::filesystem;

inline constexpr auto primes_to_rad = [](auto angle) { return angle * lab::constants<decltype(angle)>::pi / (60 * 180); };

template<typename value_type = double, typename Inputs>
void analyze(Inputs&& input_paths, stdf::path output_path, std::string_view graph_title)
{
	static_assert(std::floating_point<value_type>);
	static_assert(std::same_as<stdr::range_value_t<Inputs>, stdf::path>);

	using estimate_t = lab::estimate_t<value_type>;

	constexpr value_type zero = 0;

	estimate_t alpha;
	{
		struct
		{
			auto value_at(value_type alpha0, value_type delta_alpha) const { return alpha0 + delta_alpha; }
			auto derivative_at(value_type, value_type)               const { return stdv::repeat(value_type(1), 2); }
		} constexpr alpha_fn;
		constexpr value_type sensitivity = primes_to_rad(value_type(5) / 12);

		estimate_t alpha0(zero, lab::analyze_sample(std::array{ zero, zero, sensitivity }).mean_variance());
		estimate_t delta_alpha(primes_to_rad(value_type(45)), sensitivity * sensitivity / 12);

		alpha = lab::estimate(alpha_fn, std::array{alpha0, delta_alpha});

		std::print(
			"alpha0 = {:.5f} rad\n"
			"delta alpha = {:.5f} rad\n"
			"alpha = {:.5f} rad\n",
			alpha0,
			delta_alpha,
			alpha);
	}

	std::vector<std::pair<estimate_t /* t */, estimate_t /* v */>> regression_data;
	for (
		auto [ti_sample, previous_ti] = std::make_tuple<std::vector<value_type>, estimate_t>({}, { zero, zero });
		auto [i, input_path] : stdv::zip(stdv::iota(0uz), input_paths))
	{
		{
			std::ifstream input(input_path);
			if (!input)
				throw std::runtime_error(std::format("Cannot open {} for reading.", input_path.string()));
			input.exceptions(input.badbit);
			ti_sample.assign_range(stdv::istream<value_type>(input));
		}
		auto ti = lab::analyze_sample(ti_sample);
		
		std::print(
			"{}\n"
			"t{} sample (s): {}\n"
			"E(t): {:.5f} s\n"
			"s^2(t): {:.3} s^2\n"
			"s(t): {:.5f} s\n"
			"Skewness: {:.3f}\n"
			"Kurtosis: {:.3f}\n"
			"Size: {}\tMin: {:.4f} s\tMax: {:.4f} s\n\n",
			input_path.filename().string(),
			i, ti_sample,
			ti.mean_estimate(),
			ti.variance(),
			ti.stddev(),
			ti.skewness(),
			ti.kurtosis(),
			ti.size(), ti.min(), ti.max()
		);
		struct
		{
			auto value_at(value_type prev, value_type ti) const { return (prev + ti) / 2; }
			auto derivative_at(value_type, value_type)    const { return stdv::repeat(value_type(0.5), 2); }
		} constexpr t_fn;
		struct
		{
			auto      value_at(value_type dx, value_type prev, value_type ti) const { return dx / (ti - prev); }
			auto derivative_at(value_type dx, value_type prev, value_type ti) const
			{
				value_type
					inv_delta = 1 / (ti - prev),
					term = dx * inv_delta * inv_delta;
				return std::array{ inv_delta, term, -term };
			}
		} constexpr v_fn;
		constexpr estimate_t delta_x(0.1, 0.001 * 0.001 / 12);

		estimate_t ti_mean_estimate = ti.mean_estimate();
		regression_data.push_back(
			{
				lab::estimate(t_fn, std::array{previous_ti, ti_mean_estimate}),
				lab::estimate(v_fn, std::array{delta_x, previous_ti, ti_mean_estimate})
			});
		previous_ti = ti_mean_estimate;
	}
	std::print("t_{{i-1, i}} (s)\t\t\tv_i (ms^-1)\n");
	for (auto [x, y] : regression_data)
		std::print("{:.5f}\t:\t{:.4f}\n", x, y);

	auto regression_result = lab::regression(regression_data | stdv::transform([](auto x) {return std::pair(x.first.value(), x.second.value()); }));

	std::print(
		"\nLinear regression: v = ({:.4f}) * t + {:.4f}  (adjusted R^2: {:.5f})\n\n",
		regression_result.slope_estimate(), regression_result.intercept_estimate(), regression_result.R2()
	);

	{
		struct
		{
			auto      value_at(value_type b, value_type a) const { return b / std::sin(a); }
			auto derivative_at(value_type b, value_type a) const
			{
				value_type sin = std::sin(a);
				return std::array{ 1 / sin, -b * std::cos(a) / (sin * sin) };
			}
		} constexpr g0_fn;

		auto g0 = lab::estimate(g0_fn, std::array{regression_result.slope_estimate(), alpha});
		std::print(
			"g0 = {:.2f} ms^-2\n"
			"Compatibility with g*: {:.2f}\n\n",
			g0,
			std::abs(g0.value() - 9.806) / std::sqrt(g0.variance() + 0.001 * 0.001)
		);
	}
	{
		double
			x_min = regression_data.front().first.value(),
			x_max = regression_data.back().first.value(),
			y_min = regression_data.front().second.value(),
			y_max = regression_data.back().second.value(),
			x_margin = (x_max - x_min) / 24,
			y_margin = (y_max - y_min) / 24;
		lab::plotter2d plotter(graph_title, { "t (s)", x_min - x_margin, x_max + x_margin }, { "v (ms^{-1})", y_min - y_margin, y_max + y_margin });
		plotter.plot_regression(regression_data);
		plotter.plot_line(regression_result.slope(), regression_result.intercept(), x_min - x_margin, x_max + x_margin);
		plotter.save_as(output_path);
	}
}

int main()
{	
	stdf::path base_path = "";
	auto no_weight = stdv::iota(5, 12) | stdv::transform([&](int i) { return base_path / std::format("no_weight_40_{}0.txt", i);});
	auto weight = stdv::iota(5, 12) | stdv::transform([&](int i) { return base_path / std::format("weight_40_{}0.txt", i);});
	
	try
	{
		analyze(no_weight, base_path / "no_weight.png", "Slitta Scarica");
		std::print("\n\n\n\n");
		analyze(weight, base_path / "weight.png", "Slitta Carica");
	}
	catch (std::exception const& e)
	{
		std::print("{}", e.what());
	}
}

