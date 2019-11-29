#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/LU>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

using namespace Eigen;

constexpr int param_num = 21;
constexpr int state_dim = 6;
constexpr int input_dim = 1;
constexpr double g = 9.8;
constexpr int default_max_itr = 20;

constexpr double robotmass = 100.0;
constexpr double controldt = 0.01;

constexpr double px0 = -0.3;
constexpr double vx0 = 1.0;
constexpr double py0 = -0.1;
constexpr double vy0 = 0.3;
constexpr double pz0 = 1.0;
constexpr double vz0 = 0.1;

constexpr double px1 = 0.3;
constexpr double vx1 = 1.0;
constexpr double py1 = -0.1;
constexpr double vy1 = -0.3;
constexpr double pz1 = 1.0;
constexpr double vz1 = 0.0;

constexpr double steptime = 0.6;

constexpr double wgp = 1e4;
constexpr double wgv = 1e2;
constexpr double wp = 1e-1;
constexpr double wv = 1e-6;
constexpr double wi = 1e-6;

class DDPCoGGenerator {
public:
  void init(void)
  {
    xstart = VectorXd::Zero(state_dim);
    xgoal = VectorXd::Zero(state_dim);
    xstart << px0, vx0, py0, vy0, pz0, vz0;
    xgoal << px1, vx1, py1, vy1, pz1, vz1;
    T = steptime;
    dt = controldt;
    wgoalpos = wgp;
    wgoalvel = wgv;
    wpos = wp;
    wvel = wv;
    winput = wi;
    mass = robotmass;
    max_itr = default_max_itr;
  }

  int readConfig(std::string fname)
  {
    std::ifstream ifs(fname);
    std::string reading_line_buffer;
    if (ifs.fail()) {
      init();
      return -1;
    }
    std::vector<std::string> param_string_list;
    while (std::getline(ifs, reading_line_buffer)) {
      const char delimiter = ' ';
      std::string separated_string_buffer;
      std::istringstream line_separater(reading_line_buffer);
      std::getline(line_separater, separated_string_buffer, delimiter);
      param_string_list.push_back(separated_string_buffer);
    }
    if (param_string_list.size() < param_num) {
      init();
      return -1;
    }
    xstart = VectorXd::Zero(state_dim);
    xgoal = VectorXd::Zero(state_dim);
    T = std::stod(param_string_list[0]);
    dt = std::stod(param_string_list[1]);
    xstart(0) = std::stod(param_string_list[2]);
    xstart(1) = std::stod(param_string_list[3]);
    xstart(2) = std::stod(param_string_list[4]);
    xstart(3) = std::stod(param_string_list[5]);
    xstart(4) = std::stod(param_string_list[6]);
    xstart(5) = std::stod(param_string_list[7]);
    xgoal(0) = std::stod(param_string_list[8]);
    xgoal(1) = std::stod(param_string_list[9]);
    xgoal(2) = std::stod(param_string_list[10]);
    xgoal(3) = std::stod(param_string_list[11]);
    xgoal(4) = std::stod(param_string_list[12]);
    xgoal(5) = std::stod(param_string_list[13]);
    wgoalpos = std::stod(param_string_list[14]);
    wgoalvel = std::stod(param_string_list[15]);
    wpos = std::stod(param_string_list[16]);
    wvel = std::stod(param_string_list[17]);
    winput = std::stod(param_string_list[18]);
    mass = std::stod(param_string_list[19]);
    max_itr = std::stod(param_string_list[20]);
    std::cout << "mass: " << mass << std::endl;
    std::cout << "winput: " << winput << std::endl;
    std::cout << "xstart: ";
    for (int i = 0; i < state_dim; i++) std::cout << xstart(i) << " ";
    std::cout << std::endl;
    std::cout << "xgoal: ";
    for (int i = 0; i < state_dim; i++) std::cout << xgoal(i) << " ";
    std::cout << std::endl;
    return 0;
  }

  VectorXd calc_next(VectorXd state, VectorXd input)
  {
    double ax, ay, az;
    VectorXd tempu = input;
    if (tempu(0) < 0.0) tempu(0) = 0.0;
    if (state(4) == 0.0) {
      ax = 0.0;
      ay = 0.0;
    } else {
      ax = state(0) * tempu(0) / (state(4) * mass);
      ay = state(2) * tempu(0) / (state(4) * mass);
    }
    az = (tempu(0) - mass * g) / mass;
    VectorXd ret_state(state_dim);
    ret_state(0) = state(0) + dt * state(1) + 0.5 * dt * dt * ax;
    ret_state(1) = state(1) + dt * ax;
    ret_state(2) = state(2) + dt * state(3) + 0.5 * dt * dt * ay;
    ret_state(3) = state(3) + dt * ay;
    ret_state(4) = state(4) + dt * state(5) + 0.5 * dt * dt * az;
    ret_state(5) = state(5) + dt * az;
    return ret_state;
  }


  VectorXd ref_traj(double time)
  {
    VectorXd ret(state_dim);
    for (int i = 0; i < state_dim; i++) {
      if (i % 2 == 0) {
        ret(i) = (xgoal(i) - xstart(i)) / T * time + xstart(i);
      } else {
        ret(i) = 0.0;
      }
    }
    return ret;
  }

  double costFunction(VectorXd x, VectorXd u, double time)
  {
    MatrixXd Wx = MatrixXd::Zero(state_dim, state_dim);
    for (int i = 0; i < state_dim; i++) {
      if (i % 2 == 0) {
        if (i != 2) Wx(i, i) = wpos;
        else Wx(i, i) = /*1e-2 * */ wpos;
      } else {
        Wx(i, i) = wvel;
      }
    }
    MatrixXd Wu = MatrixXd::Zero(input_dim, input_dim);
    if (u(0) >= 0.0) {
      Wu(0, 0) = winput;
    } else {
      Wu(0, 0) = 1.0e4 * winput;
    }
    VectorXd xdiff(state_dim);
    xdiff = ref_traj(time) - x;
    return (0.5 * xdiff.transpose() * Wx * xdiff + 0.5 * u.transpose() * Wu * u)(0, 0);
    //return (0.5 * xdiff.transpose() * Wx * xdiff)(0, 0) + inputcost(u);
  }

  VectorXd lx(VectorXd x, VectorXd u, double time)
  {
    double eps = 1.0e-6;
    VectorXd ret(state_dim);
    for (int i = 0; i < state_dim; i++) {
      VectorXd tempx(state_dim);
      tempx = x;
      tempx(i) += eps;
      double val0 = costFunction(x, u, time);
      double val1 = costFunction(tempx, u, time);
      ret(i) = (val1 - val0) / eps;
    }
    return ret;
  }

  VectorXd lu(VectorXd x, VectorXd u, double time)
  {
    double eps = 1.0e-3;
    VectorXd ret(input_dim);
    for (int i = 0; i < input_dim; i++) {
      VectorXd tempu(input_dim);
      tempu = u;
      tempu(i) += eps;
      double val0 = costFunction(x, u, time);
      double val1 = costFunction(x, tempu, time);
      ret(i) = (val1 - val0) / eps;
    }
    return ret;
  }

  MatrixXd lxx(VectorXd x, VectorXd u, double time)
  {
    double eps = 1.0e-6;
    MatrixXd ret(state_dim, state_dim);
    for (int i = 0; i < state_dim; i++) {
      VectorXd tempx(state_dim);
      VectorXd val0(state_dim);
      VectorXd val1(state_dim);
      tempx = x;
      tempx(i) += eps;
      val0 = lx(x, u, time);
      val1 = lx(tempx, u, time);
      ret.block(0, i, state_dim, 1) = (val1 - val0) / eps;
    }
    return ret;
  }

  MatrixXd luu(VectorXd x, VectorXd u, double time)
  {
    double eps = 1.0e-3;
    MatrixXd ret(input_dim, input_dim);
    for (int i = 0; i < input_dim; i++) {
      VectorXd tempu(input_dim);
      VectorXd val0(input_dim);
      VectorXd val1(input_dim);
      tempu = u;
      tempu(i) += eps;
      val0 = lu(x, u, time);
      val1 = lu(x, tempu, time);
      ret.block(0, i, input_dim, 1) = (val1 - val0) / eps;
    }
    return ret;
  }

  MatrixXd lux(VectorXd x, VectorXd u, double time)
  {
    double eps = 1.0e-6;
    MatrixXd ret(input_dim, state_dim);
    for (int i = 0; i < state_dim; i++) {
      VectorXd tempx(state_dim);
      VectorXd val0(input_dim);
      VectorXd val1(input_dim);
      tempx = x;
      tempx(i) += eps;
      val0 = lu(x, u, time);
      val1 = lu(tempx, u, time);
      ret.block(0, i, input_dim, 1) = (val1 - val0) / eps;
    }
    return ret;
  }

  MatrixXd fx(VectorXd x, VectorXd u)
  {
    double eps = 1.0e-6;
    MatrixXd ret(state_dim, state_dim);
    for (int i = 0; i < state_dim; i++) {
      VectorXd tempx(state_dim);
      VectorXd val0(state_dim);
      VectorXd val1(state_dim);
      tempx = x;
      tempx(i) += eps;
      val0 = calc_next(x, u);
      val1 = calc_next(tempx, u);
      ret.block(0, i, state_dim, 1) = (val1 - val0) / eps;
    }
    return ret;
  }

  MatrixXd fu(VectorXd x, VectorXd u)
  {
    double eps = 1.0e-3;
    MatrixXd ret(state_dim, input_dim);
    for (int i = 0; i < input_dim; i++) {
      VectorXd tempu(input_dim);
      VectorXd val0(state_dim);
      VectorXd val1(state_dim);
      tempu = u;
      tempu(i) += eps;
      val0 = calc_next(x, u);
      val1 = calc_next(x, tempu);
      ret.block(0, i, state_dim, 1) = (val1 - val0) / eps;
    }
    return ret;
  }

  double goalCostFunction(VectorXd state)
  {
    VectorXd ref_goal_state(state_dim);
    ref_goal_state = ref_traj(T);
    MatrixXd W = MatrixXd::Zero(state_dim, state_dim);
    for (int i = 0; i < state_dim; i++) {
      if (i % 2 == 0) {
        W(i, i) = wgoalpos;
      } else {
        W(i, i) = wgoalvel;
      }
    }
    return (0.5 * (ref_goal_state - state).transpose() * W * (ref_goal_state - state))(0, 0);
  }

  VectorXd lastVx(VectorXd x)
  {
    double eps = 1.0e-6;
    VectorXd ret(state_dim);
    for (int i = 0; i < state_dim; i++) {
      VectorXd tempx(state_dim);
      tempx = x;
      tempx(i) += eps;
      double val0 = goalCostFunction(x);
      double val1 = goalCostFunction(tempx);
      ret(i) = (val1 - val0) / eps;
    }
    return ret;
  }

  MatrixXd lastVxx(VectorXd x)
  {
    double eps = 1.0e-6;
    MatrixXd ret(state_dim, state_dim);
    for (int i = 0; i < state_dim; i++) {
      VectorXd tempx(state_dim);
      VectorXd val0(state_dim);
      VectorXd val1(state_dim);
      tempx = x;
      tempx(i) += eps;
      val0 = lastVx(x);
      val1 = lastVx(tempx);
      ret.block(0, i, state_dim, 1) = (val1 - val0) / eps;
    }
    return ret;
  }

  void initTraj(void)
  {
    double time = 0.0;
    VectorXd x(state_dim);
    x = xstart;
    while (time < T) {
      state_list.push_back(x);
      VectorXd u(input_dim);
      u << mass * g;
      x = calc_next(x, u);
      input_list.push_back(u);
      time += dt;
    }
  }

  void calcTraj(void)
  {
    for (int cnt = 0; cnt < max_itr; cnt++) {
      backward();
      forward();
    }
  }

  void backward(void)
  {
    size_t size = state_list.size();
    VectorXd last_state(state_dim);
    last_state = calc_next(state_list[size - 1], input_list[size - 1]);
    VectorXd Vx(state_dim);
    MatrixXd Vxx(state_dim, state_dim);
    Vx = lastVx(last_state);
    Vxx = lastVxx(last_state);

    kvec_list.clear();
    Kmat_list.clear();
    for (int i = size - 1; i >= 0; i--) {
      double time = T * i / size;
      VectorXd x(state_dim);
      VectorXd u(input_dim);
      x = state_list[i];
      u = input_list[i];

      VectorXd Qx(state_dim);
      VectorXd Qu(input_dim);
      MatrixXd Qxx(state_dim, state_dim);
      MatrixXd Qux(input_dim, state_dim);
      MatrixXd Quu(input_dim, input_dim);

      VectorXd lxvec(state_dim);
      VectorXd luvec(input_dim);
      MatrixXd fxmat(state_dim, state_dim);
      MatrixXd fumat(state_dim, input_dim);
      MatrixXd lxxmat(state_dim, state_dim);
      MatrixXd luxmat(input_dim, state_dim);
      MatrixXd luumat(input_dim, input_dim);

      lxvec = lx(x, u, time);
      luvec = lu(x, u, time);
      lxxmat = lxx(x, u, time);
      luxmat = lux(x, u, time);
      luumat = luu(x, u, time);

      fxmat = fx(x, u);
      fumat = fu(x, u);

      Qx = lxvec + fxmat.transpose() * Vx;
      Qu = luvec + fumat.transpose() * Vx;
      Qxx = lxxmat + fxmat.transpose() * Vxx * fxmat;
      Qux = luxmat + fumat.transpose() * Vxx * fxmat;
      Quu = luumat + fumat.transpose() * Vxx * fumat;

      FullPivLU< MatrixXd > Quulu(Quu);

      VectorXd kvec(input_dim);
      MatrixXd Quu_inv(input_dim, input_dim);
      Quu_inv = Quulu.inverse();
      kvec = - (Quu_inv * Qu);
      MatrixXd Kmat(input_dim, state_dim);
      Kmat = - (Quu_inv * Qux);
      kvec_list.push_back(kvec);
      Kmat_list.push_back(Kmat);

      Vx = Qx - Kmat.transpose() * Quu * kvec;
      Vxx = Qxx - Kmat.transpose() * Quu * Kmat;
    }
  }

  void forward(void)
  {
    VectorXd old_state(state_dim);
    old_state = state_list[0];
    int size = state_list.size();
    for (int i = 0; i < size; i++) {
      input_list[i] = input_list[i] + kvec_list[size - 1 -i] + Kmat_list[size - 1 - i] * (state_list[i] - old_state);
      if (i == (size - 1)) break;
      old_state = state_list[i + 1];
      state_list[i + 1] = calc_next(state_list[i], input_list[i]);
    }
  }

  void output(void)
  {
    std::ofstream ofs("traj.dat");
    for (int i = 0; i < state_list.size(); i++) {
      VectorXd tempx(state_dim);
      VectorXd tempu(input_dim);
      tempx = state_list[i];
      tempu = input_list[i];
      ofs << i << " ";
      for (int j = 0; j < state_dim; j++) {
        ofs << tempx(j) << " ";
      }
      for (int j = 0; j < input_dim; j++) {
        ofs << tempu(j) << " ";
      }
      ofs << std::endl;
    }
  }

  void cogoutput(void)
  {
    std::ofstream ofs("cogtraj.l");
    ofs << "(setq *cog-traj* (list " << std::endl;
    for (int i = 0; i < state_list.size(); i++) {
      ofs << "(list ";
      VectorXd tempx(state_dim);
      VectorXd tempu(input_dim);
      tempx = state_list[i];
      tempu = input_list[i];
      ofs << i * dt << " ";
      for (int j = 0; j < state_dim; j++) {
        ofs << tempx(j) << " ";
      }
      for (int j = 0; j < input_dim; j++) {
        ofs << tempu(j) << " ";
      }
      ofs << ")" << std::endl;
    }
    ofs << "))" << std::endl;
  }

  void generate(void)
  {
    initTraj();
    calcTraj();
    output();
  }

  std::vector<VectorXd> state_list;
  std::vector<VectorXd> input_list;
  double dt;
  double T;
  double mass;
  double wgoalpos;
  double wgoalvel;
  double wpos;
  double wvel;
  double winput;
  double max_itr;
  VectorXd xstart;
  VectorXd xgoal;

  std::vector<VectorXd> kvec_list;
  std::vector<MatrixXd> Kmat_list;

};

int main(void)
{
  DDPCoGGenerator ddp;
  ddp.readConfig("param.txt");
  ddp.generate();
  return 0;
}
