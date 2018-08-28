#include <fstream>

#include <boost/unordered_map.hpp>

#include <sys/time.h>

#include <iostream>

#include <methods/ilqr.h>

#ifdef SFML_FOUND
#include <SFML/Graphics.hpp>
#endif

#include <VirtualRobot/Robot.h>
#include <VirtualRobot/RobotNodeSet.h>
#include <VirtualRobot/IK/DifferentialIK.h>

#include <VirtualRobot/XML/RobotIO.h>

using namespace OPTCTRL;
using namespace arma;



class LoadReduce
{
public:
    LoadReduce(int time, int timesteps) :
        T(time), N(timesteps)
    {
        dt = T/N;

        robot = VirtualRobot::RobotIO::loadRobot("/common/homes/students/stocker/armarx/Armar6RT/data/Armar6RT/robotmodel/Armar6-SH/Armar6-SH.xml");
        rns_left = robot->getRobotNodeSet("LeftArm");

    }

    int T, N;
    double dt;

    VirtualRobot::RobotNodeSetPtr rns_left;
    VirtualRobot::RobotPtr robot;


    vec sysfcn(const vec& x, const vec& u)
    {
        vec xnew(x.n_rows, fill::zeros);

        xnew = x + 0.1 * dt * u;

        return xnew;
    }

    double costfcn(const vec& x, const vec& u)
    {

        Eigen::VectorXf jointValues(x.size());

        for (int i = 0; i < x.size(); i++) {
            jointValues[i] = x[i];
        }

        rns_left->setJointValues(jointValues);

        VirtualRobot::RobotNodePtr tcp_left = rns_left->getTCP();

        VirtualRobot::DifferentialIKPtr ik;

        ik.reset(new VirtualRobot::DifferentialIK(rns_left, rns_left->getRobot()->getRootNode()));

        Eigen::MatrixXf jacobi;

        jacobi = ik->getJacobianMatrix(tcp_left);

        Eigen::VectorXf F_ext(6);
        F_ext << 0, 0, 1, 0, 0, 1;

        double cost = F_ext.transpose() * jacobi * jacobi.transpose() * F_ext;

        cost += 0.5 * dot(u.t(), u);

        return cost;

        //return 0;

    }

    double costfinal(const vec& x)
    {

        Eigen::Vector3f y_T(3);
        y_T << -50, 1000, 1000;

        Eigen::VectorXf jointValues(x.size());

        for (int i = 0; i < x.size(); i++) {
            jointValues[i] = x[i];
        }

        rns_left->setJointValues(jointValues);

        VirtualRobot::RobotNodePtr tcp_left = rns_left->getTCP();

        Eigen::Vector3f pos = tcp_left->getPositionInRootFrame();

        double cost = 0.5 * (y_T - pos).transpose() * (y_T - pos);

        return cost;
    }

#ifdef SFML_FOUND
    bool runSys(const vecVec& U, const vec& x0)
    {
        vec x = x0;
        double width = 1000;
        double height = 500;
        sf::RenderWindow window(sf::VideoMode(width, height), "Double Pendulum");



        double L1 = 100 * l1;
        double L2 = 100 * l2;
        sf::RectangleShape rod0(sf::Vector2f(L1+2, 5));
        rod0.setFillColor(sf::Color::Red);
        sf::RectangleShape rod1(sf::Vector2f(L2, 5));
        rod1.setFillColor(sf::Color::Red);
        sf::RectangleShape cart(sf::Vector2f(120, 50));
        cart.setFillColor(sf::Color::Green);


        int i = 0;
        while (window.isOpen())
        {
            sf::Event event;
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                    window.close();
            }

            window.clear();

            // draw ...
            double th0 = width/2 + 150 * x(0);
            double th1 = x(1) - M_PI/2;
            double th2 = x(2) - M_PI/2;

            cart.setPosition(sf::Vector2f(th0-60, height/2-25));
            rod0.setPosition(sf::Vector2f(th0, height/2));
            rod1.setPosition(sf::Vector2f(th0 + L1 * cos(th1), height/2 + L1 * sin(th1)));

//            std::cout << "th0: " << th0 << " th1: " << th1 << " "

            rod0.setRotation(th1 * 180/M_PI);
            rod1.setRotation(th2 * 180/M_PI);

            if(U.size() > i)
            {
                x = sysfcn(x, U[i]);
            }
            i++;

            window.draw(cart);
            window.draw(rod0);
            window.draw(rod1);

            window.display();
            usleep(dt * 1000000);
        }

        return true;
    }

#else
    bool runSys(const vecVec& U)
    {
        return false;
    }
#endif


    mat sfcn_x_d() {

        mat A(8, 8, fill::eye);

        return A;
    }

    mat sfcn_u_d() {

        mat A(8, 8, fill::eye);

        return dt * A;
    }

    vec cfcn_u_d(const vec& u) {

        return u;
    }

    mat cfcn_uu_d() {

        mat A(8, 8, fill::eye);

        return A;
    }


};


int main(int argc, char* argv[])
{
    LoadReduce model(10, 500);
    ilqrFcnVec sysfcn = std::bind(&LoadReduce::sysfcn, &model, std::placeholders::_1, std::placeholders::_2);
    ilqrFcnVal costfcn = std::bind(&LoadReduce::costfcn, &model, std::placeholders::_1, std::placeholders::_2);
    ilqrFcn0 costfinalfcn = std::bind(&LoadReduce::costfinal, &model, std::placeholders::_1);

    iLQR ilqr(600, 8, 8, sysfcn, costfcn, costfinalfcn);

    vec x0(8, fill::zeros);
//    x0[0] = 0;
//    x0[1] = 1;
//    x0[2] = 0;
//    x0[3] = 1;
//    x0[4] = 0;
//    x0[5] = 1;
//    x0[6] = 0;
//    x0[7] = 1;

    ilqr.eps_u = 0.01;
    ilqr.eps_x = 0.01;
    ilqr.x0 = x0;
    ilqr.convThreshold = 1e-6;

    ilqr.sfcn_x = std::bind(&LoadReduce::sfcn_x_d, &model);
    ilqr.sfcn_u = std::bind(&LoadReduce::sfcn_u_d, &model);
    ilqr.cfcn_u = std::bind(&LoadReduce::cfcn_u_d, &model, std::placeholders::_1);
    ilqr.cfcn_uu = std::bind(&LoadReduce::cfcn_uu_d, &model);

    vecVec U = ilqr.run(100);

    //model.runSys(U, x0);


}


