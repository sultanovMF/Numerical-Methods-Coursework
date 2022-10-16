#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iomanip>
#include <tuple>
const double pi = 3.1415926535897931;
// Библиотеки для вывода графиков/табличек
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>
#include <implot/implot.h> 

class Mesh {
	// qnnlwnl
public:
	Mesh(const double L, const double T, const double dx, const double dt) {
		len_x = floor(L / dx) + 1;
		len_t = floor(T / dt) + 1;

		mesh.resize(len_x * len_t);
		for (int j = 0; j < len_t; ++j) {
			for (int i = 0; i < len_x; ++i) {
				auto& [x, t, u] = mesh[i + j * len_x];
				x = i * dx;
				t = j * dt;
			}
		}

	}

	std::tuple<double, double, double>& operator()(const int i, const int j) {
		// TODO проверять вылезаем ли за сетку
		return mesh[i + j * len_x];
	}
	std::tuple<double, double, double> operator()(const int i, const int j) const {
		// TODO проверять вылезаем ли за сетку
		return mesh[i + j * len_x];
	}

	unsigned int getLenX() const {
		return len_x;
	}

	unsigned int getLenT() const {
		return len_t;
	}
private:
	std::vector<std::tuple<double, double, double>> mesh; // where tuple contains x, t, u values
	unsigned int len_x;
	unsigned int len_t;
};

void solve_wave(
	Mesh& mesh,
	const double C,
	const std::function<double(double)> f_0,  
	const std::function<double(double)> f_l,   
	const std::function<double(double)> g_1,  
	const std::function<double(double)> g_2) {
	/// Solving wave equation u_tt = a^2 u_xx
	//* C : Courant numebr
	//* u(0, t) = f_0 
	//* u(l, t) = f_l
	//* u(x,0)	= g_1
	//* D[2](u)(x,0) = g_2
	const unsigned int len_x = mesh.getLenX();
	const unsigned int len_t = mesh.getLenT();

	const double dx = std::get<0>(mesh(1, 0)) - std::get<0>(mesh(0, 0));
	const double dt = std::get<1>(mesh(0, 1)) - std::get<1>(mesh(0, 0));


	// Setting initial condition
	for (int i = 0; i < len_x; ++i) {
		auto& [x, t, u] = mesh(i, 0);
		u = g_1(x);
	}

	for (int i = 0; i < len_x; ++i) {
		auto& [x, t, u] = mesh(i, 1);
		// TODO можно улучшить если использовать разложение g_2 в ряд тейлора
		u = g_1(x) + dt * g_2(x);
	}

	{
		auto& [x, t, u] = mesh(0, 1);
		u = f_0(t);
	}

	{
		auto& [x, t, u] = mesh(len_x - 1, 1);
		u = f_l(t);
	}

	const double C2 = C * C;
	for (int j = 1; j < len_t - 1; ++j) {
		{
			// setting boundary condition
			auto& [x_left, t_left, u_left] = mesh(0, j + 1);
			auto& [x_right,t_right, u_right] = mesh(len_x - 1, j + 1);
			u_left = f_0(t_left);
			u_right = f_l(t_right);
		}

		for (int i = 1; i < len_x - 1; ++i) {
			double& u_next			= std::get<2>(mesh(i, j + 1));
			const double u_left		= std::get<2>(mesh(i - 1, j));
			const double u_right	= std::get<2>(mesh(i + 1, j));
			const double u_center	= std::get<2>(mesh(i, j));
			const double u_bottom	= std::get<2>(mesh(i, j - 1));
			
			u_next = C2 * (u_right - 2 * u_center + u_left) + 2 * u_center - u_bottom;
		}
	}
}


double zero(double) {
	return 0.;
}

double foo(double x) {
	if (0 <= x <= 1) return x;
	else if (1 <= x <= 2) return 2-x;
}

int main() {
	double L = 2;
	double T = 2;
	double dx = 0.005;
	double dt = 0.001;
	double a = 1;
	double C = dt / dx * a;

	Mesh mesh(L, T, dx, dt);
	solve_wave(mesh, C, zero, zero, zero, foo);
	
	std::vector<double> numerical_x;
	std::vector<double> numerical_u;

	const int j = 3;
	int plot_point = 10;
	const int ratio = round(mesh.getLenX() / (plot_point));
	plot_point++;
	for (int i = 0; i < plot_point; ++i) {
		auto [x, t, u] = mesh(i * ratio, j);
		numerical_x.push_back(x);
		numerical_u.push_back(u);
	}


	std::vector<double> target_x;
	std::vector<double> target_u;

	for (int i = 0; i < plot_point; ++i) {
		auto [x, t, u] = mesh(i * ratio, j);
		target_x.push_back(x);

		double result = 0;
		for (int n = 1; n < 10000; n+=2) {
			result += 8. / (pi * pi * n * n) * sin(pi * n * t / 2) * sin(pi * n * x / 2);
			result -= 8. / (pi * pi * (n+1) * (n+1)) * sin(pi * (n + 1) * t / 2) * sin(pi * (n + 1) * x / 2);
		}

		target_u.push_back(result);
	}


    sf::RenderWindow window(sf::VideoMode(1280, 720), "ImGui + SFML = <3");
    window.setFramerateLimit(60);
    ImGui::SFML::Init(window);

    sf::Clock deltaClock;

    ImGui::CreateContext();
    ImPlot::CreateContext();


    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            ImGui::SFML::ProcessEvent(event);

            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        ImGui::SFML::Update(window, deltaClock.restart());

        static bool use_work_area = true;
        static ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
        // ImGui::ShowDemoWindow();
        const ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(use_work_area ? viewport->WorkPos : viewport->Pos);
        ImGui::SetNextWindowSize(use_work_area ? viewport->WorkSize : viewport->Size);



        ImGui::Begin("Coursework");
		if (ImPlot::BeginPlot("Approx")) {
			ImPlot::PlotLine("Approx", &numerical_x[0], &numerical_u[0], plot_point);
			
			ImPlot::EndPlot();
		}
		if (ImPlot::BeginPlot("Target")) {
			ImPlot::PlotLine("Target", &target_x[0], &target_u[0], plot_point);
			ImPlot::EndPlot();
		}
        ImGui::End();

        window.clear();
        ImGui::SFML::Render(window);
        window.display();
    }


    ImGui::SFML::Shutdown();
	return 0;
}