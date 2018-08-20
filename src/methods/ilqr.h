#ifndef ILQR_H
#define ILQR_H

#include <armadillo>
#include <functional>
using namespace arma;
using namespace std::placeholders;

namespace OPTCTRL
{
    typedef std::vector<mat> vecMat;
    typedef std::vector<vec> vecVec;

    typedef std::function<double(const vec&) > ilqrFcn0;
    typedef std::function<vec(const vec&) > ilqrFcn0Vec;
    typedef std::function<mat(const vec&) > ilqrFcn0Mat;

    typedef std::function<double(const vec&, const vec&) > ilqrFcnVal;
    typedef std::function<vec(const vec&, const vec&) > ilqrFcnVec;
    typedef std::function<mat(const vec&, const vec&) > ilqrFcnMat;

    class iLQR
    {
    public:
        iLQR(){}
        iLQR(int timeSteps, int stateDim, int ctrlDim,
             ilqrFcnVec sysfcn, ilqrFcnVal costfcn, ilqrFcn0 costfinalfcn)
        {
            this->stateDim = stateDim;
            this->ctrlDim = ctrlDim;
            T = timeSteps;
            sfcn = sysfcn;
            cfcn = costfcn;
            cffcn = costfinalfcn;

            eps_x = 1e-2;
            eps_u = 1e-2;

            sfcn_x = std::bind(&iLQR::sfcn_x_d, this, _1, _2);
            sfcn_u = std::bind(&iLQR::sfcn_u_d, this, _1, _2);

            cfcn_x = std::bind(&iLQR::cfcn_x_d, this, _1, _2);
            cfcn_u = std::bind(&iLQR::cfcn_u_d, this, _1, _2);
            cfcn_xx = std::bind(&iLQR::cfcn_xx_d, this, _1, _2);
            cfcn_ux = std::bind(&iLQR::cfcn_ux_d, this, _1, _2);
            cfcn_uu = std::bind(&iLQR::cfcn_uu_d, this, _1, _2);

            cffcn_x = std::bind(&iLQR::cffcn_x_d, this, _1);
            cffcn_xx = std::bind(&iLQR::cffcn_xx_d, this, _1);

            lamb = 1.0;
            lamb_factor = 10;
            lamb_max = 1000;
            convThreshold = 1e-3;

            x0 = vec(stateDim, fill::zeros);
        }

        ilqrFcnVec sfcn;
        ilqrFcnVal cfcn;
        ilqrFcn0 cffcn;
        ilqrFcn0Vec cffcn_x;
        ilqrFcn0Mat cffcn_xx;

        int T;

        // partial derivations
        ilqrFcnMat sfcn_x;
        ilqrFcnMat sfcn_u;


        ilqrFcnVec cfcn_x;
        ilqrFcnVec cfcn_u;
        ilqrFcnMat cfcn_xx;
        ilqrFcnMat cfcn_ux;
        ilqrFcnMat cfcn_uu;

        vecMat fx;
        vecMat fu;
        vecMat lxx;
        vecMat lux;
        vecMat luu;
        vecVec lu;
        vecVec lx;

        vecVec k;
        vecMat K;

        double V;
        vec Vx, Vx0;
        mat Vxx, Vxx0;

        vecVec U;
        vecVec X;
        int stateDim;
        int ctrlDim;
        double eps_x;
        double eps_u;
        double lamb;
        double lamb_factor;
        double lamb_max;
        double convThreshold;
        vec x0;

        double rollout()
        {
            if(U.size() < T-1)
            {
                std::cerr << "The control signal seq is smaller than the required time steps, appending random control value" << std::endl;

                while (U.size() < T-1)
                {
                    vec u(ctrlDim, fill::zeros);
                    U.push_back(u);
                }
            }

            vec x = x0;
            double cost = 0;
            fx.clear();
            fu.clear();
            lx.clear();
            lu.clear();
            lxx.clear();
            lux.clear();
            luu.clear();
            k.clear();
            K.clear();

            for(size_t i = 0; i < U.size(); ++i)
            {
                vec u = U[i];
                fx.push_back(sfcn_x(x, u));
                fu.push_back(sfcn_u(x, u));
                lx.push_back(cfcn_x(x, u));
                lu.push_back(cfcn_u(x, u));
                lxx.push_back(cfcn_xx(x, u));
                lux.push_back(cfcn_ux(x, u));
                luu.push_back(cfcn_uu(x, u));

                k.push_back(vec(ctrlDim,fill::zeros));
                K.push_back(mat(ctrlDim, stateDim, fill::zeros));

                cost += cfcn(x, u);

                if(X.size() == T)
                {
                    x = X[i+1];
                }
                else
                {
                    x = sfcn(x, u);
                    X.push_back(x);
                }
            }

            double c = cffcn(x);
            V = c;
            Vx0 = cffcn_x(x);
            Vxx0 = cffcn_xx(x);
            cost += c;

            return cost;
        }

        vecVec Unew;
        vecVec Xnew;
        double backpass()
        {
            vec Qx, Qu;
            mat Qxx, Qux, Quu;

            Vx = Vx0;
            Vxx = Vxx0;
            for(int t = T-2; t > 0; t--)
            {
                Qx = lx[t] + fx[t].t() * Vx;
                Qu = lu[t] + fu[t].t() * Vx;

                Qxx = lxx[t] + fx[t].t() * Vxx * fx[t];
                Qux = lux[t] + fu[t].t() * Vxx * fx[t];
                Quu = luu[t] + fu[t].t() * Vxx * fu[t];

                vec eigval;
                mat eigvec;
                eig_sym(eigval, eigvec, Quu);

                for(int i = 0; i < eigval.n_rows; ++i)
                {
                    if(eigval(i) < 0)
                    {
                        eigval(i) = 0;
                    }
                }
                eigval = eigval + lamb;
                eigval = 1 / eigval;
                mat emat = diagmat(eigval);
                mat Quu_i = eigvec * emat * eigvec.t();

                k[t] = - Quu_i * Qu;
                K[t] = - Quu_i * Qux;

                Vx = Qx - K[t].t() * Quu * k[t];
                Vxx = Qxx - K[t].t() * Quu * K[t];
            }


            vec xnew = x0;
            Unew = U;
            Xnew = X;
            double newcost = 0;
            for(int t = 0; t < T-1; t++)
            {
                Unew[t] = U[t] + k[t] + K[t] * (xnew - X[t]);
                newcost += cfcn(xnew, Unew[t]);
                xnew = sfcn(xnew, Unew[t]);
                Xnew.push_back(xnew);
            }
            newcost += cffcn(xnew);


            return newcost;

        }



        vecVec run(vec& xres, int maxEpoch=100)
        {
            bool createNewTraj = true;
            double oldcost;
            int i = 0;
            for(; i < maxEpoch; ++i)
            {
                if(createNewTraj)
                {
                    oldcost = rollout();
                    createNewTraj = false;
                }

                double newcost = backpass();

                if(newcost < oldcost)
                {
                    lamb /= lamb_factor;
                    U = Unew;
                    X = Xnew;
                    createNewTraj = true;
                    std::cout << "Iter: " << i << "\t Cost: " << newcost << std::endl;
                    if(i > 0 && (fabs(newcost - oldcost) / newcost) < convThreshold)
                    {
                        std::cout << "converged at Iter: " << i << "\t Cost: " << newcost << std::endl;
                        break;
                    }
                }
                else
                {
                    lamb *= lamb_factor;
                    if(lamb > lamb_max)
                    {
                        std::cout << "lambda is greater than the maximal value at Iter: " << i << "\t Cost: " << newcost << std::endl;
                        break;
                    }
                }

            }

            std::cout << "[Final Report] Iter: " << i << "\t Cost: " << oldcost << std::endl;
            xres = X[X.size()-1];
            return U;
        }

    private:
        mat sfcn_x_d(const vec& x, const vec& u)
        {
            mat A(x.n_rows, x.n_rows, fill::zeros);
            for(int i = 0; i < x.n_rows; ++i)
            {
                vec x_inc = x;
                x_inc[i] += eps_x;
                x_inc = sfcn(x_inc, u);
                vec x_dec = x;
                x_dec[i] -= eps_x;
                x_dec = sfcn(x_dec, u);

                A.col(i) = (x_inc - x_dec) / (2 * eps_x);
            }

            return A;
        }

        mat sfcn_u_d(const vec& x, const vec& u)
        {
            mat A(x.n_rows, u.n_rows, fill::zeros);
            for(int i = 0; i < u.n_rows; ++i)
            {
                vec u_inc = u;
                u_inc[i] += eps_u;
                vec x_inc = sfcn(x, u_inc);
                vec u_dec = u;
                u_dec[i] -= eps_u;
                vec x_dec = sfcn(x, u_dec);
                A.col(i) = (x_inc - x_dec) / (2 * eps_u);
            }

            return A;
        }

        vec cffcn_x_d(const vec& x)
        {
            vec A(x.n_rows,fill::zeros);
            for(int i = 0; i < x.n_rows; ++i)
            {
                vec x_inc = x;
                x_inc[i] += eps_x;
                double c_inc = cffcn(x_inc);
                vec x_dec = x;
                x_dec[i] -= eps_x;
                double c_dec = cffcn(x_dec);
                A(i) = (c_inc - c_dec) / (2 * eps_x);
            }
            return A;
        }


        mat cffcn_xx_d(const vec& x)
        {
            mat A(x.n_rows, x.n_rows, fill::zeros);
            for(int i = 0; i < x.n_rows; ++i)
            {
                for(int j = 0; j < x.n_rows; ++j)
                {
                    vec x_inc_inc = x;
                    x_inc_inc[i] += eps_x;
                    x_inc_inc[j] += eps_x;
                    double c_inc_inc = cffcn(x_inc_inc);
                    vec x_inc_dec = x;
                    x_inc_dec[i] += eps_x;
                    x_inc_dec[j] -= eps_x;
                    double c_inc_dec = cffcn(x_inc_dec);
                    vec x_dec_inc = x;
                    x_dec_inc[i] -= eps_x;
                    x_dec_inc[j] += eps_x;
                    double c_dec_inc = cffcn(x_dec_inc);
                    vec x_dec_dec = x;
                    x_dec_dec[i] -= eps_x;
                    x_dec_dec[j] -= eps_x;
                    double c_dec_dec = cffcn(x_dec_dec);

                    A(i,j) = (c_inc_inc - c_inc_dec - c_dec_inc + c_dec_dec) / 4 * eps_x * eps_x;
                }
            }

            return A;
        }


        vec cfcn_x_d(const vec& x, const vec& u)
        {
            vec A(x.n_rows,fill::zeros);
            for(int i = 0; i < x.n_rows; ++i)
            {
                vec x_inc = x;
                x_inc[i] += eps_x;
                double c_inc = cfcn(x_inc, u);
                vec x_dec = x;
                x_dec[i] -= eps_x;
                double c_dec = cfcn(x_dec, u);
                A(i) = (c_inc - c_dec) / (2 * eps_x);
            }
            return A;
        }     

        vec cfcn_u_d(const vec& x, const vec& u)
        {
            vec A(u.n_rows,fill::zeros);
            for(int i = 0; i < u.n_rows; ++i)
            {
                vec u_inc = u;
                u_inc[i] += eps_u;
                double c_inc = cfcn(x, u_inc);
                vec u_dec = u;
                u_dec[i] -= eps_u;
                double c_dec = cfcn(x, u_dec);
                A(i) = (c_inc - c_dec) / (2 * eps_u);
            }
            return A;
        }

        mat cfcn_xx_d(const vec& x, const vec& u)
        {
            mat A(x.n_rows, x.n_rows, fill::zeros);
            for(int i = 0; i < x.n_rows; ++i)
            {
                for(int j = 0; j < x.n_rows; ++j)
                {
                    vec x_inc_inc = x;
                    x_inc_inc[i] += eps_x;
                    x_inc_inc[j] += eps_x;
                    double c_inc_inc = cfcn(x_inc_inc, u);
                    vec x_inc_dec = x;
                    x_inc_dec[i] += eps_x;
                    x_inc_dec[j] -= eps_x;
                    double c_inc_dec = cfcn(x_inc_dec, u);
                    vec x_dec_inc = x;
                    x_dec_inc[i] -= eps_x;
                    x_dec_inc[j] += eps_x;
                    double c_dec_inc = cfcn(x_dec_inc, u);
                    vec x_dec_dec = x;
                    x_dec_dec[i] -= eps_x;
                    x_dec_dec[j] -= eps_x;
                    double c_dec_dec = cfcn(x_dec_dec, u);

                    A(i,j) = (c_inc_inc - c_inc_dec - c_dec_inc + c_dec_dec) / 4 * eps_x * eps_x;
                }
            }

            return A;
        }

        mat cfcn_ux_d(const vec& x, const vec& u)
        {
            mat A(u.n_rows, x.n_rows, fill::zeros);
            for(int i = 0; i < u.n_rows; ++i)
            {
                for(int j = 0; j < x.n_rows; ++j)
                {
                    vec u_inc = u;
                    u_inc[i] += eps_u;
                    vec x_inc = x;
                    x_inc[j] += eps_x;
                    vec u_dec = u;
                    u_dec[i] -= eps_u;
                    vec x_dec = x;
                    x_dec[j] -= eps_x;
                    double c_inc_inc = cfcn(x_inc, u_inc);
                    double c_inc_dec = cfcn(x_dec, u_inc);
                    double c_dec_inc = cfcn(x_inc, u_dec);
                    double c_dec_dec = cfcn(x_dec, u_dec);

                    A(i,j) = (c_inc_inc - c_inc_dec - c_dec_inc + c_dec_dec) / 4 * eps_x * eps_u;
                }
            }
            return A;
        }

        mat cfcn_uu_d(const vec& x, const vec& u)
        {
            mat A(u.n_rows, u.n_rows, fill::zeros);
            for(int i = 0; i < u.n_rows; ++i)
            {
                for(int j = 0; j < u.n_rows; ++j)
                {
                    vec u_inc_inc = u;
                    u_inc_inc[i] += eps_u;
                    u_inc_inc[j] += eps_u;
                    double c_inc_inc = cfcn(x, u_inc_inc);
                    vec u_inc_dec = u;
                    u_inc_dec[i] += eps_u;
                    u_inc_dec[j] -= eps_u;
                    double c_inc_dec = cfcn(x, u_inc_dec);
                    vec u_dec_inc = u;
                    u_dec_inc[i] -= eps_u;
                    u_dec_inc[j] += eps_u;
                    double c_dec_inc = cfcn(x, u_dec_inc);
                    vec u_dec_dec = u;
                    u_dec_dec[i] -= eps_u;
                    u_dec_dec[j] -= eps_u;
                    double c_dec_dec = cfcn(x, u_dec_dec);

                    A(i,j) = (c_inc_inc - c_inc_dec - c_dec_inc + c_dec_dec) / 4 * eps_u * eps_u;
                }
            }

            return A;
        }

    };
}


#endif
