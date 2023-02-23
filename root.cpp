export module lab:root;

import :core;

import <TAxis.h>;
import <TCanvas.h>;
import <TGraphErrors.h>;
import <TF1.h>;
import <TH1F.h>;

export namespace lab
{
	struct plotter2d
	{
	private:
		TCanvas _canvas;
		std::stack<std::unique_ptr<TGraph>> _graphs;
	public:

		struct axis_params
		{
			std::string_view label;
			double min, max;
		};

		plotter2d(std::string_view title, axis_params x, axis_params y)
		{
			_canvas.SetName("");
			_canvas.SetCanvasSize(4096, 2160);
			auto frame = _canvas.DrawFrame(x.min, y.min, x.max, y.max, title.data());

			frame->GetXaxis()->SetTitle(x.label.data());
			frame->GetYaxis()->SetTitle(y.label.data());
		}

		template<typename Range>
		void plot_regression(Range&& data)
		{
			static_assert(stdr::range<Range>);
			std::vector<double> x, y;
			std::vector<double> ex, ey;
			if constexpr (stdr::sized_range<Range>)
			{
				size_t size = stdr::size(data);
				x.reserve(size);
				y.reserve(size);
				ex.reserve(size);
				ey.reserve(size);
			}
			for (auto [first, second] : data)
			{
				x.push_back(first.value());
				y.push_back(second.value());
				ex.push_back(first.stddev());
				ey.push_back(second.stddev());
			}
			_canvas.cd();
			_graphs.emplace(std::make_unique<TGraphErrors>(static_cast<Int_t>(x.size()), x.data(), y.data(), ex.data(), ey.data()))
				->Draw("P");
		}

		void plot_line(double m, double q, double x_min, double x_max)
		{
			std::array<double, 2> x = { x_min, x_max }, y = { m * x_min + q, m * x_max + q };
			_canvas.cd();
			_graphs.emplace(std::make_unique<TGraph>(2, x.data(), y.data()))
				->Draw("C");
		}

		void save_as(stdf::path const& path)
		{
			_canvas.SaveAs(path.string().c_str());
		}
	};
}
