#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/LU>
#include <vector>
#include <fstream>

using namespace Eigen;

double start_xpos = -0.2;
double start_xvel = 1.0;
double start_zpos = 0.7;
double start_zvel = 0.1;
double goal_xpos = 0.2;
double goal_xvel = 1.0;
double goal_zpos = 0.7;
double goal_zvel = -0.1;

double mass = 100.0;
double g = 9.8;

double T = 0.5;
double dt = 0.01;

const int state_num = 4;
const int input_num = 2;


VectorXd ref_traj(double time)
{
  double midz_height = start_zpos + 0.05;
  VectorXd ret(state_num);
  ret(1) = 0.0;
  ret(3) = 0.0;
  ret(0) = (goal_xpos - start_xpos) / T * time + start_xpos;
  if (time < 0.5 * T) {
    ret(2) = (midz_height - start_zpos) / (0.5 * T) * time + start_zpos;
  } else {
    ret(2) = (goal_zpos - midz_height) / (0.5 * T) * (time - 0.5 * T) + midz_height;
  }
  return ret;
}

VectorXd calc_next(VectorXd state, VectorXd input)
{
  double ax, az;
  ax = state(0) * input(0) / (state(2) * mass);
  //ax = input(1) / mass;
  az = (input(0) - mass * g) / mass;
  VectorXd ret_state(state_num);
  ret_state << (state(0) + dt * state(1) + 0.5 * dt * dt * ax),
            (state(1) + dt * ax),
            (state(2) + dt * state(3) + 0.5 * dt * dt * az),
            (state(3) + dt * az);
  return ret_state;
}

double inputcost(VectorXd u)
{
  double u0 = u(0) - mass * g;
  double w1 = 1.0e2;
  double w2 = 1.0e-6;
  double w3 = 1.0e2;
  /*
  if (u0 > 0.0) return 0.5 * u0 * u0 * w1 + w3;
  else if (u(0) < 0.0) return 0.5 * u(0) * u(0) * w2;
  else return w3 * (u(0) / (mass * g)) * (u(0) / (mass * g));
  */
  if (u0 > 0.0) {
    return 0.5 * u(0) * u(0) * 1.0e-8;
  } else if (u(0) < 0.0) {
    return 0.5 * u(0) * u(0) * 1.0e-4;
  } else {
    return 0.5 * u(0) * u(0) * 1.0e-8;
  }
}

double costFunction(VectorXd state, VectorXd u, double time)
{
  MatrixXd Wstate = MatrixXd::Zero(state_num, state_num);
  Wstate(0, 0) = 1e1;
  Wstate(1, 1) = 1e-5;
  Wstate(2, 2) = 1e1;
  Wstate(3, 3) = 1e-5;
  MatrixXd Winput = MatrixXd::Zero(input_num, input_num);
  Winput(0, 0) = 1e-6;
  Winput(1, 1) = 1e-6;
  VectorXd diff_state(state_num);
  diff_state = ref_traj(time) - state;
  return (0.5 * diff_state.transpose() * Wstate * diff_state + 0.5 * u.transpose() * Winput * u)(0, 0);
  //return (0.5 * diff_state.transpose() * Wstate * diff_state)(0, 0) + inputcost(u);
}

VectorXd lx(VectorXd state, VectorXd u, double time)
{
  double eps = 1.0e-6;
  VectorXd ret(state_num);
  for (int i = 0; i < state_num; i++) {
    VectorXd tempstate(state_num);
    tempstate = state;
    tempstate(i) += eps;
    double val0 = costFunction(state, u, time);
    double val1 = costFunction(tempstate, u, time);
    ret(i) = (val1 - val0) / eps;
  }
  return ret;
}

VectorXd lu(VectorXd state, VectorXd u, double time)
{
  double eps = 1.0e-3;
  VectorXd ret(input_num);
  for (int i = 0; i < input_num; i++) {
    VectorXd tempu(input_num);
    tempu = u;
    tempu(i) += eps;
    double val0 = costFunction(state, u, time);
    double val1 = costFunction(state, tempu, time);
    ret(i) = (val1 - val0) / eps;
  }
  return ret;
}

MatrixXd lxx(VectorXd state, VectorXd u, double time)
{
  double eps = 1.0e-6;
  MatrixXd ret(state_num, state_num);
  for (int i = 0; i < state_num; i++) {
    VectorXd tempstate(state_num);
    VectorXd tempret1(state_num);
    VectorXd tempret0(state_num);
    tempstate = state;
    tempstate(i) += eps;
    tempret0 = lx(state, u, time);
    tempret1 = lx(tempstate, u, time);
    for (int j = 0; j < state_num; j++) {
      ret(j, i) = (tempret1(j) - tempret0(j)) / eps;
    }
  }
  return ret;
}

MatrixXd luu(VectorXd state, VectorXd u, double time)
{
  double eps = 1.0e-3;
  MatrixXd ret(input_num, input_num);
  for (int i = 0; i < input_num; i++) {
    VectorXd tempu(input_num);
    VectorXd tempret1(input_num);
    VectorXd tempret0(input_num);
    tempu = u;
    tempu(i) += eps;
    tempret0 = lu(state, u, time);
    tempret1 = lu(state, tempu, time);
    for (int j = 0; j < input_num; j++) {
      ret(j, i) = (tempret1(j) - tempret0(j)) / eps;
    }
  }
  return ret;
}

MatrixXd lux(VectorXd state, VectorXd u, double time)
{
  double eps = 1.0e-6;
  MatrixXd ret(input_num, state_num);
  for (int i = 0; i < state_num; i++) {
    VectorXd tempstate(state_num);
    VectorXd tempret0(input_num);
    VectorXd tempret1(input_num);
    tempstate = state;
    tempstate(i) += eps;
    tempret0 = lu(state, u, time);
    tempret1 = lu(tempstate, u, time);
    for (int j = 0; j < input_num; j++) {
      ret(j, i) = (tempret1(j) - tempret0(j)) / eps;
    }
  }
  return ret;
}

MatrixXd fx(VectorXd state, VectorXd u)
{
  double eps = 1.0e-6;
  MatrixXd ret(state_num, state_num);
  for (int i = 0; i < state_num; i++) {
    VectorXd tempstate(state_num);
    VectorXd tempret0(state_num);
    VectorXd tempret1(state_num);
    tempstate = state;
    tempstate(i) += eps;
    tempret0 = calc_next(state, u);
    tempret1 = calc_next(tempstate, u);
    for (int j = 0; j < state_num; j++) {
      ret(j, i) = (tempret1(j) - tempret0(j)) / eps;
    }
  }
  return ret;
}

MatrixXd fu(VectorXd state, VectorXd u)
{
  double eps = 1.0e-3;
  MatrixXd ret(state_num, input_num);
  for (int i = 0; i < input_num; i++) {
    VectorXd tempu(input_num);
    VectorXd tempret0(state_num);
    VectorXd tempret1(state_num);
    tempu = u;
    tempu(i) += eps;
    tempret0 = calc_next(state, u);
    tempret1 = calc_next(state, tempu);
    for (int j = 0; j < state_num; j++) {
      ret(j, i) = (tempret1(j) - tempret0(j)) / eps;
    }
  }
  return ret;
}

double lastCostFunction(VectorXd state)
{
  VectorXd ref_goal_state(state_num);
  ref_goal_state = ref_traj(T);
  MatrixXd W(state_num, state_num);
  W = Matrix<double, state_num, state_num>::Zero();
  W(0, 0) = 1e4;
  W(1, 1) = 10.0;
  W(2, 2) = 1e4;
  W(3, 3) = 10.0;

  return (0.5 * (ref_goal_state - state).transpose() * W * (ref_goal_state - state))(0, 0);
}

VectorXd lastVx(VectorXd state)
{
  double eps = 1.0e-6;
  VectorXd ret(state_num);
  for (int i = 0; i < state_num; i++) {
    VectorXd tempstate(state_num);
    tempstate = state;
    tempstate(i) += eps;
    double val0 = lastCostFunction(state);
    double val1 = lastCostFunction(tempstate);
    ret(i) = (val1 - val0) / eps;
  }
  return ret;
}

MatrixXd lastVxx(VectorXd state)
{
  double eps = 1.0e-6;
  MatrixXd ret(state_num, state_num);
  for (int i = 0; i < state_num; i++) {
    VectorXd tempstate(state_num);
    VectorXd tempret0(state_num);
    VectorXd tempret1(state_num);
    tempstate = state;
    tempstate(i) += eps;
    tempret0 = lastVx(state);
    tempret1 = lastVx(tempstate);
    for (int j = 0; j < state_num; j++) {
      ret(j, i) = (tempret1(j) - tempret0(j)) / eps;
    }
  }
  return ret;
}

void ddptest(void)
{
  std::vector<VectorXd> state_list;
  std::vector<VectorXd> input_list;


  {
    VectorXd x(state_num);
    x << start_xpos, start_xvel, start_zpos, start_zvel;
    double time = 0.0;
    while (time < T) {
      state_list.push_back(x);
      VectorXd u(input_num);
      u << mass * g, 0.0;
      x = calc_next(x, u);
      input_list.push_back(u);
      time += dt;
    }
  }

  int size = state_list.size();

  std::ofstream ofs1("first.dat");
  for (int i = 0; i < size; i++) {
    VectorXd temp(state_num);
    temp = state_list[i];
    ofs1 << i << " " << temp(0) << " " << temp(2) << " " << input_list[i](0) << " " << input_list[i](1) << std::endl;
  }

  for (int cnt = 0; cnt < 50; cnt++) {
    std::cout << "cnt: " << cnt << std::endl;
    VectorXd last_state(state_num);
    last_state = calc_next(state_list[size - 1], input_list[size - 1]);
    VectorXd Vx(state_num);
    MatrixXd Vxx(state_num, state_num);
    Vx = lastVx(last_state);
    Vxx = lastVxx(last_state);

    std::vector<VectorXd> kvec_list;
    std::vector<MatrixXd> Kmat_list;
    for (int i = size - 1; i >= 0; i--) {
      double time = T * i / size;
      VectorXd x(state_num);
      VectorXd u(input_num);
      x = state_list[i];
      u = input_list[i];

      VectorXd Qx(state_num);
      VectorXd Qu(input_num);
      MatrixXd Qxx(state_num, state_num);
      MatrixXd Qux(input_num, state_num);
      MatrixXd Quu(input_num, input_num);

      VectorXd lxvec(state_num);
      VectorXd luvec(input_num);
      MatrixXd fxmat(state_num, state_num);
      MatrixXd fumat(state_num, input_num);
      MatrixXd lxxmat(state_num, state_num);
      MatrixXd luxmat(input_num, state_num);
      MatrixXd luumat(input_num, input_num);

      lxvec = lx(x, u, time);
      luvec = lu(x, u, time);
      lxxmat = lxx(x, u ,time);
      luxmat = lux(x, u, time);

      /*
      std::cout << "x: ";
      for (int ii = 0; ii < state_num; ii++) std::cout << x(ii) << " ";
      std::cout << std::endl;
      std::cout << "u: ";
      for (int ii = 0; ii < input_num; ii++) std::cout << u(ii) << " ";
      std::cout << std::endl;
      */
      luumat = luu(x, u, time);
      //std::cout << "luu_mat: ";
      //std::cout << std::endl;
      //for (int ii = 0; ii < input_num; ii++) {
      //  for (int jj = 0; jj < input_num; jj++) {
      //    std::cout << luumat(ii, jj) << " ";
      //  }
      //  std::cout << std::endl;
      //}
      //std::cout << std::endl;

      fxmat = fx(x, u);
      fumat = fu(x, u);

      Qx = lxvec + fxmat.transpose() * Vx;
      Qu = luvec + fumat.transpose() * Vx;
      Qxx = lxxmat + fxmat.transpose() * Vxx * fxmat;
      Qux = luxmat + fumat.transpose() * Vxx * fxmat;
      Quu = luumat + fumat.transpose() * Vxx * fumat;

      FullPivLU< MatrixXd > Quulu(Quu);
      //std::cout << "rank: " << lu.rank() << " " << input_num << std::endl;

      VectorXd kvec(input_num);
      MatrixXd Quu_inv(2, 2);
      //Quu_inv = Quu.inverse();
      Quu_inv = Quulu.inverse();
      kvec = - (Quu_inv * Qu);
      MatrixXd Kmat(input_num, state_num);
      Kmat = - (Quu_inv * Qux);
      kvec_list.push_back(kvec);
      Kmat_list.push_back(Kmat);
      
      /*
      for (int ii = 0; ii < state_num; ii++) std::cout << Vx(ii) << " ";
      std::cout << std::endl;
      std::cout << std::endl;
      for (int ii = 0; ii < state_num; ii++) {
        for (int jj = 0; jj < state_num; jj++) {
          std::cout << Vxx(ii, jj) << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      */
      //std::cout << diff_state(0) << " " << diff_state(2) << " " << Kmat(0, 0) << std::endl;
      //std::cout << "Qx: ";
      //for (int ii = 0; ii < state_num; ii++) std::cout << Qx(ii) << " ";
      //std::cout << std::endl;
      //std::cout << "lxvec: ";
      //for (int ii = 0; ii < state_num; ii++) std::cout << lxvec(ii) << " ";
      //std::cout << std::endl;
      //std::cout << "luvec: ";
      //for (int ii = 0; ii < input_num; ii++) std::cout << luvec(ii) << " ";
      //std::cout << std::endl;
      //std::cout << "Vx: ";
      //for (int ii = 0; ii < state_num; ii++) std::cout << Vx(ii) << " ";
      //std::cout << std::endl;
      //std::cout << "kvec: ";
      //for (int ii = 0; ii < input_num; ii++) std::cout << kvec(ii) << " ";
      //std::cout << std::endl;
      //std::cout << "Kmat: ";
      //std::cout << std::endl;
      //for (int ii = 0; ii < input_num; ii++) {
      //  for (int jj = 0; jj < state_num; jj++) {
      //    std::cout << Kmat(ii, jj) << " ";
      //  }
      //  std::cout << std::endl;
      //}
      //std::cout << std::endl;
      //std::cout << "Quu: ";
      //std::cout << std::endl;
      //for (int ii = 0; ii < input_num; ii++) {
      //  for (int jj = 0; jj < input_num; jj++) {
      //    std::cout << Quu(ii, jj) << " ";
      //  }
      //  std::cout << std::endl;
      //}
      //std::cout << std::endl;
      //std::cout << "fumat: ";
      //std::cout << std::endl;
      //for (int ii = 0; ii < state_num; ii++) {
      //  for (int jj = 0; jj < input_num; jj++) {
      //    std::cout << fumat(ii, jj) << " ";
      //  }
      //  std::cout << std::endl;
      //}
      //std::cout << std::endl;
      //std::cout << "Vxx: ";
      //std::cout << std::endl;
      //for (int ii = 0; ii < state_num; ii++) {
      //  for (int jj = 0; jj < state_num; jj++) {
      //    std::cout << Vxx(ii, jj) << " ";
      //  }
      //  std::cout << std::endl;
      //}
      //std::cout << std::endl;
      //std::cout << std::endl;

      Vx = Qx - Kmat.transpose() * Quu * kvec;
      Vxx = Qxx - Kmat.transpose() * Quu * Kmat;
    }

    VectorXd old_state(state_num);
    old_state = state_list[0];
    for (int i = 0; i < size; i++) {
      input_list[i] = input_list[i] + kvec_list[size - 1 - i] + Kmat_list[size - 1 - i] * (state_list[i] - old_state);
      if (i == size - 1) break;
      old_state = state_list[i + 1];
      state_list[i + 1] = calc_next(state_list[i], input_list[i]);
    }
  }

  std::ofstream ofs2("last.dat");
  for (int i = 0; i < size; i++) {
    VectorXd temp(state_num);
    temp = state_list[i];
    ofs2 << i << " " << temp(0) << " " << temp(1) << " " << temp(2) << " " << temp(3) << " " << input_list[i](0) << " " << input_list[i](1) << std::endl;
  }

  std::ofstream ofs3("cogtraj.l");
  ofs3 << "(setq *cog-traj* (list " << std::endl;
  for (int i = 0; i < size; i++) {
    VectorXd temp(state_num);
    temp = state_list[i];
    ofs3 << "(list ";
    ofs3 << (i * dt / T) << " " << temp(0) << " " << temp(2) << ") " << std::endl;
  }
  ofs3 << "))";
}

void eigentest(void)
{
  VectorXd x(state_num);
  VectorXd u(input_num);
  x << 0.69, 1.0, 1.0, 0.0;
  u << 392.1, 0.0;

  VectorXd luvec(input_num);
  MatrixXd lumat(input_num, input_num);
  luvec = lu(x, u, 0);
  lumat = luu(x, u, 0);
  std::cout << "luvec: ";
  for (int i = 0; i < input_num; i++) {
    std::cout << luvec(i) << " ";
  }
  std::cout << std::endl;
  std::cout << "lumat: ";
  for (int i = 0; i < input_num; i++) {
    for (int j = 0; j < input_num; j++) {
      std::cout << lumat(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  x = ref_traj(T);
  for (int i = 0; i < state_num; i++) std::cout << x(i) << " ";
  std::cout << std::endl;

}

int main(void)
{
//  eigentest();
  ddptest();
  return 0;
}
