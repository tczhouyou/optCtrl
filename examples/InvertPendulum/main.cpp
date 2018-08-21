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

class invPendulum
{
public:
    invPendulum(double mass, double length, double gravity, double friction, double timestep = 1e-3) : m(mass), l(length), g(gravity), mu(friction), dt(timestep){}

    double m, l, g, mu, dt;

    vec sysfcn(const vec& x, const vec& u)
    {
        double dx1 = x[1];
        double dx2 = g * sin(x[0]) / l  - mu * x[1] / (m * l * l) + u[0] / (m * l * l);

        vec xnew = x;
        xnew[0] += dx1 * dt;
        xnew[1] += dx2 * dt;

        return xnew;
    }

    double costfcn(const vec& x, const vec& u)
    {
        return 0.5 * u[0] * u[0] * 1e-5 + 0.5 * (x[0] * x[0] + x[1] * x[1]);
    }

    double costfinal(const vec& x)
    {
        return 0.5 * (x[0] * x[0] + x[1] * x[1]);
    }

#ifdef SFML_FOUND
    bool runSys(const vecVec& U, const vec& x0)
    {
        vec x = x0;
        sf::RenderWindow window(sf::VideoMode(500, 500), "Inverted Pendulum");
        sf::RectangleShape line(sf::Vector2f(100, 5));
        line.setFillColor(sf::Color::Red);
        line.setPosition(250, 250);

        int i = 0;
        float theta0 = x[0] - M_PI/2;
        while (window.isOpen())
        {
            sf::Event event;
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                    window.close();
            }

            window.clear();

            float theta = x[0] - M_PI/2;
            theta = theta * 180 / M_PI;
            line.setRotation(theta);

            if(U.size() > i)
            {
                x = sysfcn(x, U[i]);
            }
            theta0 = theta;
            i++;


            window.draw(line);
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
    invPendulum model(1, 1, 9.8, 0.01, 1e-2);
    ilqrFcnVec sysfcn = std::bind(&invPendulum::sysfcn, &model, _1, _2);
    ilqrFcnVal costfcn = std::bind(&invPendulum::costfcn, &model, _1, _2);
    ilqrFcn0 costfinalfcn = std::bind(&invPendulum::costfinal, &model, _1);

    iLQR ilqr(200, 2, 1, sysfcn, costfcn, costfinalfcn);

    vec x0(2, fill::zeros);
    x0[0] = M_PI;

    ilqr.eps_u = 0.01;
    ilqr.eps_x = 0.01;
    ilqr.x0 = x0;
    ilqr.convThreshold = 1e-6;

    vecVec Xres;
    vecVec U = ilqr.mpc(Xres, 1000, 200);

    model.runSys(U, x0);


}


