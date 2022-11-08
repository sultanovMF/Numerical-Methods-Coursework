//#include <iostream>
//#include <cmath>
//#include <functional>
//#include <algorithm>
//#include <vector>
//#include <cmath>
//#include <iomanip>
//#include <tuple>
//const double pi = 3.1415926535897931;
//// Библиотеки для вывода графиков/табличек
//#include <SFML/Graphics/RenderWindow.hpp>
//#include <SFML/System/Clock.hpp>
//#include <SFML/Window/Event.hpp>
//#include <imgui-SFML.h>
//#include <imgui.h>
//#include <implot/implot.h> 
//
//class Mesh {
//public:
//	Mesh(const double L, const double T, const double dx, const double dt) {
//		len_x = floor(L / dx) + 1;
//		len_t = floor(T / dt) + 1;
//
//		mesh.resize(len_x * len_t);
//		for (int j = 0; j < len_t; ++j) {
//			for (int i = 0; i < len_x; ++i) {
//				auto& [x, t, u] = mesh[i + j * len_x];
//				x = i * dx;
//				t = j * dt;
//			}
//		}
//
//	}
//
//	std::tuple<double, double, double>& operator()(const int i, const int j) {
//		// TODO проверять вылезаем ли за сетку
//		return mesh[i + j * len_x];
//	}
//	std::tuple<double, double, double> operator()(const int i, const int j) const {
//		// TODO проверять вылезаем ли за сетку
//		return mesh[i + j * len_x];
//	}
//
//	unsigned int getLenX() const {
//		return len_x;
//	}
//
//	unsigned int getLenT() const {
//		return len_t;
//	}
//private:
//	std::vector<std::tuple<double, double, double>> mesh; // where tuple contains x, t, u values
//	unsigned int len_x;
//	unsigned int len_t;
//
//	int a;
//};
//
//void solve_wave(
//	Mesh& mesh,
//	const double C,
//	const std::function<double(double)> f_0,  
//	const std::function<double(double)> f_l,   
//	const std::function<double(double)> g_1,  
//	const std::function<double(double)> g_2) {
//	/// Solving wave equation u_tt = a^2 u_xx
//	//* C : Courant numebr
//	//* u(0, t) = f_0 
//	//* u(l, t) = f_l
//	//* u(x,0)	= g_1
//	//* D[2](u)(x,0) = g_2
//	const unsigned int len_x = mesh.getLenX();
//	const unsigned int len_t = mesh.getLenT();
//
//	const double dx = std::get<0>(mesh(1, 0)) - std::get<0>(mesh(0, 0));
//	const double dt = std::get<1>(mesh(0, 1)) - std::get<1>(mesh(0, 0));
//
//	// Setting initial condition
//	for (int i = 0; i < len_x; ++i) {
//		auto& [x, t, u] = mesh(i, 0);
//		u = g_1(x);
//	}
//
//	for (int i = 0; i < len_x; ++i) {
//		auto& [x, t, u] = mesh(i, 1);
//		// TODO можно улучшить если использовать разложение g_2 в ряд тейлора
//		u = g_1(x) + dt * g_2(x);
//	}
//	// setting boundary condition for first layer
//	{
//		auto& [x, t, u] = mesh(0, 1);
//		u = f_0(t);
//	}
//
//	{
//		auto& [x, t, u] = mesh(len_x - 1, 1);
//		u = f_l(t);
//	}
//
//	const double C2 = C * C;
//	for (int j = 1; j < len_t - 1; ++j) {
//		{
//			// setting boundary condition
//			auto& [x_left, t_left, u_left] = mesh(0, j + 1);
//			auto& [x_right,t_right, u_right] = mesh(len_x - 1, j + 1);
//			u_left = f_0(t_left);
//			u_right = f_l(t_right);
//		}
//
//		for (int i = 1; i < len_x - 1; ++i) {
//			double& u_next			= std::get<2>(mesh(i, j + 1));
//			const double u_left		= std::get<2>(mesh(i - 1, j));
//			const double u_right	= std::get<2>(mesh(i + 1, j));
//			const double u_center	= std::get<2>(mesh(i, j));
//			const double u_bottom	= std::get<2>(mesh(i, j - 1));
//			
//			u_next = C2 * (u_right - 2 * u_center + u_left) + 2 * u_center - u_bottom;
//		}
//	}
//}
//
//
//class Mesh3D {
//public:
//	Mesh3D(const double L_x, const double L_y, const double L_z, const double L_t, 
//		const double dl, const double dt) {
//		
//		len_x = ceil(L_x / dl) + 1;
//		len_y = ceil(L_y / dl) + 1;
//		len_z = ceil(L_z / dl) + 1;
//		len_t = ceil(L_t / dt) + 1;
//
//		mesh.resize(len_x * len_y * len_z * len_t);
//
//		for (int i = 0; i < len_t; ++i) {
//			for (int j = 0; j < len_x; ++j) {
//				for (int k = 0; k < len_y; ++k) {
//					for (int n = 0; n < len_z; n++) {
//						auto& [x, y, z, t, u] = 
//							mesh[i + j * len_x + k * len_x * len_y + n * len_x * len_y * len_z];
//						x = j * dl;
//						y = k * dl;
//						z = n * dl;
//						t = i * dt;
//					}
//				}
//			}
//		}
//	}
//
//	unsigned int getLenX() const {
//		return len_x;
//	}
//	unsigned int getLenY() const {
//		return len_y;
//	}
//	unsigned int getLenZ() const {
//		return len_z;
//	}
//	unsigned int getLenT() const {
//		return len_t;
//	}
//
//
//	std::tuple<double, double, double, double, double>& operator()
//		(const int i, const int j, const int k, const int n) {
//		// TODO проверять вылезаем ли за сетку
//		return mesh[i + j * len_x + k * len_x * len_y + n * len_x * len_y * len_z];
//	}
//	std::tuple<double, double, double, double, double> operator()
//		(const int i, const int j, const int k, const int n) const {
//		// TODO проверять вылезаем ли за сетку
//		return mesh[i + j * len_x + k * len_x * len_y + n * len_x * len_y * len_z];
//	}
//private:
//	std::vector<std::tuple<double, double, double, double, double>> mesh; // where tuple contains x, y, z, t, u values
//	unsigned int len_x;
//	unsigned int len_y;
//	unsigned int len_z;
//	unsigned int len_t;
//};
//
//
//void solve_wave_3d(
//	Mesh3D& mesh,
//	const double C,
//	const std::function<double(double)> f_0,
//	const std::function<double(double)> f_l,
//	const std::function<double(double, double, double)> g_1,
//	const std::function<double(double, double, double)> g_2) {
//
//	const unsigned int len_x = mesh.getLenX();
//	const unsigned int len_y = mesh.getLenY();
//	const unsigned int len_z = mesh.getLenZ();
//	const unsigned int len_t = mesh.getLenT();
//
//	const double dl = std::get<0>(mesh(1, 0, 0, 0)) - std::get<0>(mesh(0, 0, 0, 0));
//	const double dt = std::get<1>(mesh(0, 0, 0, 1)) - std::get<1>(mesh(0, 0, 0, 1));
//
//	
// Build mesh (with i j k) support
// Setting initial  condition (x, y, z)
// Setting boundary condition (t)
// Evaluateing (only prev u)
// 
// 
// 
//	// Setting initial condition
//	{
//		for (int i = 0; i < len_x; ++i) {
//			for (int j = 0; j < len_y; ++j) {
//				for (int k = 0; k < len_z; ++k) {
//					auto& [x, y, z, t, u] = mesh(i, j, k, 0);
//					u = g_1(x, y, z);
//				}
//			}
//		}
//		for (int i = 1; i < len_x; ++i) {
//			for (int j = 1; j < len_y; ++j) {
//				for (int k = 1; k < len_z; ++k) {
//					auto& [x, y, z, t, u] = mesh(i, j, k, 1);
//					// TODO можно улучить если увеличить класс гладкости
//					u = g_1(x, y, z) + dt * g_2(x, y, z);
//				}
//			}
//		}
//	}
//	//* Set boundary condition
//	auto set_all_boundaries = [&](double time) {
//		auto set_boundary = [&mesh, f_0, time](double len_x, double len_y, double len_z,
//			bool take_x, bool take_y, bool take_z) {
//				for (int i = 0 - !take_x; i < len_x * take_x; ++i * take_x) {
//					for (int j = 0 - !take_y; j < len_y * take_y; ++j * take_y) {
//						for (int k = 0 - !take_z; k < len_z * take_z; ++k * take_z) {
//							// if (i == -1 && j == -1; && k == -1;)
//							auto& [x, y, z, t, u] = mesh(
//								i * take_x + !take_x * len_x,
//								j * take_y + !take_y * len_y,
//								k * take_z + !take_z * len_z, time);
//							u = f_0(t);
//						}
//					}
//				}
//		};
//
//		set_boundary(0, 0, 0, true, true, false);
//		set_boundary(0, 0, 0, true, false, true);
//		set_boundary(0, 0, 0, false, true, false);
//		set_boundary(0, 0, len_z - 1, true, true, false);
//		set_boundary(0, len_y - 1, 0, true, false, true);
//		set_boundary(len_x - 1, 0, 0, false, true, false);
//	};
//
//	set_all_boundaries(1);
//	
//	const double C2 = C * C;
//	for (int n = 1; n < len_t - 1; ++n) {
//		set_all_boundaries(n+1);
//		for (int i = 1; i < len_x - 1; ++i) {
//			for (int j = 1; j < len_y - 1; ++j) {
//				for (int k = 1; k < len_z - 1; ++k) {
//					double& u_next = std::get<4>(mesh(i, j, k , n+1));
//					{
//						const double u_left = std::get<4>(mesh(i - 1, j, k, n));
//						const double u_right = std::get<4>(mesh(i + 1, j, k, n));
//						const double u_center = std::get<4>(mesh(i, j, k, n));
//						u_next += C2 * (u_right - 2 * u_center + u_left);
//					}
//					{
//						const double u_left = std::get<4>(mesh(i, j - 1, k, n));
//						const double u_right = std::get<4>(mesh(i, j + 1, k, n));
//						const double u_center = std::get<4>(mesh(i, j, k, n));
//						u_next += C2 * (u_right - 2 * u_center + u_left);
//					}
//					{
//						const double u_left = std::get<4>(mesh(i, j, k - 1, n));
//						const double u_right = std::get<4>(mesh(i, j, k + 1, n));
//						const double u_center = std::get<4>(mesh(i, j, k, n));
//						u_next += C2 * (u_right - 2 * u_center + u_left);
//					}
//					{
//						const double u_center = std::get<4>(mesh(i, j, k, n));
//						const double u_bottom = std::get<4>(mesh(i, j, k, n - 1));
//						u_next += 2 * u_center - u_bottom;
//					}
//				}
//			}
//		}
//	}
//}
//
//double zero(double) {
//	return 0.;
//}
//
//int main() {
//	double L = 1;
//	double T = 1;
//	double dx = 0.005;
//	double dt = 0.001;
//	double a = 0.333;
//	double C = dt / dx * a;
//
//	Mesh mesh(L, T, dx, dt);
//	solve_wave(mesh, C, zero, zero, [](double x) { return sin(pi * x); }, zero);
//	
//	std::vector<double> numerical_x;
//	std::vector<double> numerical_u;
//
//	const int j = 10;
//	int plot_point = 100;
//	const int ratio = round(mesh.getLenX() / (plot_point));
//	//plot_point++;
//	for (int i = 0; i < plot_point; ++i) {
//		auto [x, t, u] = mesh(i * ratio, j);
//		numerical_x.push_back(x);
//		numerical_u.push_back(u);
//	}
//
//
//	std::vector<double> target_x;
//	std::vector<double> target_u;
//
//	for (int i = 0; i < plot_point; ++i) {
//		auto [x, t, u] = mesh(i * ratio, j);
//		target_x.push_back(x);
//
//		double result = 0;
//		//for (int n = 1; n < 10000; n+=2) {
//		//	/*result += 8. / (pi * pi * n * n) * sin(pi * n * t / 2) * sin(pi * n * x / 2);
//		//	result -= 8. / (pi * pi * (n+1) * (n+1)) * sin(pi * (n + 1) * t / 2) * sin(pi * (n + 1) * x / 2);*/
//		//}
//		result = cos(pi * t) * sin(pi * x);
//		target_u.push_back(result);
//	}
//
//	//Mesh3D mesh(L, L, L, T, dx, dt);
//
//	//std::vector <double> target_x;
//	//std::vector <double> target_u;
//	//int plot_point = 100;
//	//const int ratio = round(mesh.getLenX() / (plot_point));
//	//plot_point++;
//	//// зафиксируем y, z, t
//	//for (int i = 0; i < plot_point; ++i) {
//	//	auto [x, y, z, t, u] = mesh(i * ratio, 1, 1, 1);
//	//	target_x.push_back(x);
//	//	double result = 0;
//	//	for (int l = 0; l < 100; ++l) {
//	//		for (int m = 0; m < 100; ++m) {
//	//			for (int n = 0; n < 100; ++n) {
//	//				result += sin(l * pi * x) * sin(m * pi * y) * sin(n * pi * z) * cos(sqrt(l * l + m * m + n * n));
//	//			}
//	//		}
//	//	}
//	//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//    sf::RenderWindow window(sf::VideoMode(1280, 720), "ImGui + SFML = <3");
//    window.setFramerateLimit(60);
//    ImGui::SFML::Init(window);
//
//    sf::Clock deltaClock;
//
//    ImGui::CreateContext();
//    ImPlot::CreateContext();
//
//
//    while (window.isOpen()) {
//        sf::Event event;
//        while (window.pollEvent(event)) {
//            ImGui::SFML::ProcessEvent(event);
//
//            if (event.type == sf::Event::Closed) {
//                window.close();
//            }
//        }
//
//        ImGui::SFML::Update(window, deltaClock.restart());
//
//        static bool use_work_area = true;
//        static ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
//        // ImGui::ShowDemoWindow();
//        const ImGuiViewport* viewport = ImGui::GetMainViewport();
//        ImGui::SetNextWindowPos(use_work_area ? viewport->WorkPos : viewport->Pos);
//        ImGui::SetNextWindowSize(use_work_area ? viewport->WorkSize : viewport->Size);
//
//
//
//        ImGui::Begin("Coursework");
//		if (ImPlot::BeginPlot("Approx")) {
//			ImPlot::PlotLine("Approx", &numerical_x[0], &numerical_u[0], plot_point);
//			
//			ImPlot::EndPlot();
//		}
//		if (ImPlot::BeginPlot("Target")) {
//			ImPlot::PlotLine("Target", &target_x[0], &target_u[0], plot_point);
//			ImPlot::EndPlot();
//		}
//        ImGui::End();
//
//        window.clear();
//        ImGui::SFML::Render(window);
//        window.display();
//    }
//
//
//    ImGui::SFML::Shutdown();
//	return 0;
//}
