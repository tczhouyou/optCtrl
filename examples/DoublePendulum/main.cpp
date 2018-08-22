#include <fstream>

#include <boost/unordered_map.hpp>

#include <sys/time.h>

#include <iostream>

#include <methods/ilqr.h>

#ifdef SFML_FOUND
#include <SFML/Graphics.hpp>
#endif

using namespace OPTCTRL;
using namespace arma;



class DoublePendulumCart
{
public:
    DoublePendulumCart(double cMass, double p0Mass, double p1Mass, double l1, double l2, double g, double timestep=1e-2) :
        m0(cMass), m1(p0Mass), m2(p1Mass), l1(l1), l2(l2), g(g), dt(timestep)
    {
        d1 = m0 + m1 + m2;
        d2 = (m1/2 + m2) * l1;
        d3 = m2 * l2/2;
        d4 = (m1/3 + m2) * l1 * l1;
        d5 = m2 * l1 * l2 / 2;
        d6 = m2 * l2 * l2 / 3;
        f1 = (m1/2 + m2) * l1 * g;
        f2 = m2 * l2 * g / 2;
    }

    double m0, m1, m2, l1, l2, g, dt;
    double d1, d2, d3, d4, d5, d6, f1, f2;

    vec sysfcn(const vec& x, const vec& u)
    {
        double th1 = x(1);
        double th2 = x(2);
        double dth0 = x(3);
        double dth1 = x(4);
        double dth2 = x(5);

        vec dth(3, fill::zeros);
        dth(0) = dth0;
        dth(1) = dth1;
        dth(2) = dth2;

        mat D(3, 3, fill::zeros);
        D(0,0) = d1;
        D(1,1) = d2;
        D(2,2) = d3;
        D(0,1) = d2 * cos(th1);
        D(1,0) = D(0,1);
        D(0,2) = d3 * cos(th2);
        D(2,0) = D(0,2);
        D(1,2) = d5 * cos(th1 - th2);
        D(2,1) = D(1,2);

        mat C(3,3, fill::zeros);
        C(0,1) = -d2 * sin(th1) * dth1;
        C(0,2) = -d3 * sin(th2) * dth2;
        C(1,2) = d5 * sin(th1 - th2) * dth2;
        C(2,1) = -d5 * sin(th1 - th2) * dth1;

        vec G(3, fill::zeros);
        G(1) = -f1 * sin(th1);
        G(2) = -f2 * sin(th2);

        vec H(3, fill::zeros);
        H(0) = 1;

        vec ddth = D.i() * (H * u(0) - C * dth - G);

        vec xnew(x.n_rows, fill::zeros);
        xnew(0) = x(0) + dt * x(3);
        xnew(1) = x(1) + dt * x(4);
        xnew(2) = x(2) + dt * x(5);
        xnew(3) = x(3) + dt * ddth(0);
        xnew(4) = x(4) + dt * ddth(1);
        xnew(5) = x(5) + dt * ddth(2);

        return xnew;
    }

    double costfcn(const vec& x, const vec& u)
    {
        return u[0] * u[0] * 1e-50 + (x(1) * x(1) + x(2) * x(2));
    }

    double costfinal(const vec& x)
    {
        return (x(0) * x(0) + x(1) * x(1) + x(2) * x(2) + 1e-10 * x(5) * x(5));
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

};


int main(int argc, char* argv[])
{
    DoublePendulumCart model(1, 1, 1, 1, 1, 9.8, 1e-2);
    ilqrFcnVec sysfcn = std::bind(&DoublePendulumCart::sysfcn, &model, _1, _2);
    ilqrFcnVal costfcn = std::bind(&DoublePendulumCart::costfcn, &model, _1, _2);
    ilqrFcn0 costfinalfcn = std::bind(&DoublePendulumCart::costfinal, &model, _1);

    iLQR ilqr(300, 6, 1, sysfcn, costfcn, costfinalfcn);

    vec x0(6, fill::zeros);
    x0[0] = 0;
    x0[1] = M_PI;
    x0[2] = M_PI;

    ilqr.eps_u = 0.01;
    ilqr.eps_x = 0.01;
    ilqr.x0 = x0;
    ilqr.convThreshold = 1e-6;

    vecVec Xres;
    vecVec U = ilqr.mpc(Xres, 500, 200);

#ifdef SFML_FOUND
    model.runSys(U, x0);
#endif
    ilqr.save(U, "DoublePendulum_U");

}


