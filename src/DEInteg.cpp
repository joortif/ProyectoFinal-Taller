#include <valarray>
#include <iostream>
#include "../include/DEInteg.h"
#include "../include/SAT_Const.h"
#include "../include/Sign_.h"

using namespace std;

struct DE_STATE {
    int DE_INIT = 1;      //Restart integration
    int DE_DONE = 2;      //Successful step
    int DE_BADACC = 3;    //Accuracy requirement could not be achieved
    int DE_NUMSTEPS = 4;  //Permitted number of steps exceeded
    int DE_STIFF = 5;     //Stiff problem suspected
    int DE_INVPARAM = 6;  //Invalid input parameters
};

Matrix DEInteg(Matrix (*func)(double, Matrix), double t, double tout, double relerr, double abserr, double n_eqn, Matrix y) {
    // maxnum = 500;
    double twou  = 2*eps;
    double fouru= 4*eps;

    DE_STATE DE_STATE;
    int State_ = DE_STATE.DE_INIT;
    bool PermitTOUT = true;         // Allow integration past tout by default
    int told = 0;

    // Powers of two (two(n)=2^n)
    Matrix two(1,14);
    two(0,0)=1.0;
    two(0,1)=2.0;
    two(0,2)=4.0;
    two(0,3)=8.0;
    two(0,4)=16.0;
    two(0,5)=32.0;
    two(0,6)=64.0;
    two(0,7)=128.0;
    two(0,8)=256.0;
    two(0,9)=512.0;
    two(0,10)=1024.0;
    two(0,11)=2048.0;
    two(0,12)=4096.0;
    two(0,13)=8192.0;

    Matrix gstr(1,14);
    gstr(0,0)=1.0;
    gstr(0,1)=0.5;
    gstr(0,2)=0.0833;
    gstr(0,3)=0.0417;
    gstr(0,4)=0.0264;
    gstr(0,5)=0.0188;
    gstr(0,6)=0.0143;
    gstr(0,7)=0.0114;
    gstr(0,8)=0.00936;
    gstr(0,9)=0.00789;
    gstr(0,10)=0.00679;
    gstr(0,11)=0.00592;
    gstr(0,12)=0.00524;
    gstr(0,13)=0.00468;

    Matrix yy(n_eqn,1);
    Matrix wt(n_eqn,1);
    Matrix p(n_eqn,1);
    Matrix yp(n_eqn,1);
    Matrix phi(n_eqn,17);
    Matrix g(14,1);
    Matrix sig(14,1);
    Matrix rho(14,1);
    Matrix w(13,1);
    Matrix alpha(13,1);
    Matrix beta(13,1);
    Matrix v(13,1);
    Matrix psi_(13,1);

    // while(true)

    // Return, if output time equals input time

    if(t==tout){   // No integration
        return y;
    }

    //Test for improper parameters

    double epsilon = max(relerr,abserr);

    if (relerr < 0.0                             //Negative relative error bound
            || abserr < 0.0                      //Negative absolute error bound
            || epsilon <= 0.0                    //Both error bounds are non-positive
            || State_ > DE_STATE.DE_INVPARAM     //Invalid status flag
            || (State_ != DE_STATE.DE_INIT
            && t != told)
            ) {
        State_ = DE_STATE.DE_INVPARAM;  // Set error code
        return y;                       // Exit
    }

    // On each call set interval of integration and counter for
    // number of steps. Adjust input error tolerances to define
    // weight vector for subroutine STEP.

    double del    = tout - t;
    double absdel = fabs(del);

    double tend   = t + 100.0*del;
    if (!PermitTOUT) {
        tend = tout;
    }


    int nostep = 0;
    int kle4   = 0;
    bool stiff  = false;
    double releps = relerr/epsilon;
    double abseps = abserr/epsilon;


    bool start;
    double x;
    double delsgn;
    double h;
    bool OldPermit;

    if((State_==DE_STATE.DE_INIT) || (!OldPermit) || (delsgn*del<=0.0)){
        start = true;
        x = t;
        yy = y;
        delsgn = sign_(1.0, del);
        h = sign_(max(fouru * fabs(x), fabs(tout - x)), tout - x);
    }


    int cont=0;
    int cont2=0;
    double kold;
    double erkm1;
    double absh;
    double hnew;

    while (true) {
        if (fabs(x - t) >= absdel) {
            Matrix yout(n_eqn, 1);
            Matrix ypout(n_eqn, 1);
            g(1, 0) = 1.0;
            rho(1, 0) = 1.0;
            double hi = tout - x;
            double ki = kold + 1;

            // Initialize w[*] for computing g[*]

            // Compute g[*]
            double term = 0.0;
            for (int j = 1; j < ki; j++) {
                double psijm1 = psi_(j, 0);
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;
                for (int i = 0; i < ki + 1 - j; i++) {
                    w(i + 1, 0) = gamma * w(i + 1, 0) - eta * w(i + 2, 0);
                }
                g(j + 1, 0) = w(1, 0);
                rho(j + 1, 0) = gamma * rho(j, 0);
                term = psijm1;
            }

            // Interpolate for the solution yout and for
            // the derivative of the solution ypout

            for (int j = 0; j < ki; j++) {
                int i = ki + 1 - j - 2;
                Matrix phi_Intermediario(phi.fils(), 1);
                phi_Intermediario.print();
                for (int u = 0; u<phi.fils(); u++) {
                    phi_Intermediario(u, 0) = phi(u, i + 1);
                }
                yout = yout + phi_Intermediario * g(i + 1, 0);
                ypout = ypout + phi_Intermediario * rho(i + 1, 0);
                yout.print();
                ypout.print();
            }
            yout = y + yout * hi;
            y = yout;
            State_ = DE_STATE.DE_DONE; // Set return code
            t = tout;             // Set independent variable
            told = t;                // Store independent variable
            OldPermit = PermitTOUT;
            return y;                       // Normal exit
        }
        // If cannot go past output point and sufficiently close,
        // extrapolate and return
        if (!PermitTOUT && (fabs(tout - x) < fouru * fabs(x))) {
            h = tout - x;
            yp = func(x, yy);          // Compute derivative yp(x)
            y = yy + yp * h;                // Extrapolate vector from x to tout
            State_ = DE_STATE.DE_DONE; // Set return code
            t = tout;             // Set independent variable
            told = t;                // Store independent variable
            bool OldPermit = PermitTOUT;
            return y;                       // Normal exit
        }

        /*
         If cannot go past output point and sufficiently close,
         extrapolate and return
        if ( ~PermitTOUT && ( abs(tout-x) < fouru*abs(x) ) )
            h = tout - x;
        yp = func(x,yy);          % Compute derivative yp(x)
        y = yy + h*yp;                % Extrapolate vector from x to tout
        State_    = DE_STATE.DE_DONE; % Set return code
        t         = tout;             % Set independent variable
        told      = t;                % Store independent variable
        OldPermit = PermitTOUT;
        return;                       % Normal exit
        end

         Test for too much work
           if (nostep >= maxnum)
               State_ = DE_STATE.DE_NUMSTEPS; // Too many steps
                                                          %       if (stiff)
                   State_ = DE_STATE.DE_STIFF;% Stiffness suspected
                                                           %       end
                                                           %       y         = yy;                % Copy last step
                                                                                                              %       t         = x;
               told      = t;
               OldPermit = true;
              return;        // Weak failure exit
           end
        */
        //Limit step size, set weight vector and take a step
        h = sign_(min(fabs(h), fabs(tend - x)), h);
        for (int l = 0; l < n_eqn; l++) {
            wt(l, 0) = releps * fabs(yy(l, 0)) + abseps;
        }
        //wt.print();

        /*
         * %   Step
        %
        % Begin block 0
        %
        % Check if step size or error tolerance is too small for machine
        % precision.  If first step, initialize phi array and estimate a
        % starting step size. If step size is too small, determine an
        % acceptable one.
        %
         */

        if (fabs(h) < fouru * fabs(x)) {
            h = sign_(fouru * fabs(x), h);
            bool crash = true;
            return y;           // Exit
        }

        double p5eps = 0.5 * epsilon;
        bool crash = false;
        g(1, 0) = 1.0;
        g(2, 0) = 0.5;
        sig(1, 0) = 1.0;

        int ifail = 0;

        // If error tolerance is too small, increase it to an
        // acceptable value.

        double round = 0.0;
        for (int l = 0; l < n_eqn; l++) {
            //cout<<"----------------------"<<endl;
            //cout<<setprecision(40)<<"round: "<<round<<endl;
            //cout<<setprecision(40)<<"y(l, 0): "<<y(l, 0)<<endl;
            //cout<<setprecision(40)<<"(y(l, 0) * y(l, 0)): "<<(y(l, 0) * y(l, 0))<<endl;
            //cout<<setprecision(40)<<"wt(l, 0): "<<wt(l, 0)<<endl;
            //cout<<setprecision(40)<<"(wt(l, 0) * wt(l, 0)): "<<(wt(l, 0) * wt(l, 0))<<endl;
            //cout<<setprecision(40)<<"Division: "<<(y(l, 0) * y(l, 0)) / (wt(l, 0) * wt(l, 0))<<endl;

            round = round + (y(l, 0) * y(l, 0)) / (wt(l, 0) * wt(l, 0));
            //cout<<setprecision(40)<<"round2: "<<round<<endl;
            //cout<<"----------------------"<<endl;
        }
        round = twou * sqrt(round);
        if (p5eps < round) {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            return y;
        }
        int k;
        double hold;
        bool nornd;
        if (start) {
            // Initialize. Compute appropriate step size for first step.
            yp = func(x, y);
            yp.print();
            double sum = 0.0;
            for (int l = 0; l < n_eqn; l++) {
                phi(l, 1) = yp(l, 0);
                phi(l, 2) = 0.0;
                sum = sum + (yp(l, 0) * yp(l, 0)) / (wt(l, 0) * wt(l, 0));
            }
            sum = sqrt(sum);
            absh = fabs(h);
            if (epsilon < 16.0 * sum * h * h) {
                absh = 0.25 * sqrt(epsilon / sum);
            }
            h = sign_(max(absh, fouru * fabs(x)), h);
            hold = 0.0;
            hnew = 0.0;
            k = 1;
            kold = 0;
            start = false;
            bool phase1 = true;
            nornd = true;
            if (p5eps <= 100.0 * round) {
                nornd = false;
                for (int l = 0; l < n_eqn; l++) {
                    phi(l, 15) = 0.0;
                }
            }
        }

        /*
         * End block 0
         */

        /*
         * Repeat blocks 1, 2 (and 3) until step is successful
         */

        int ns;
        int kp1;
        int kp2;
        double km1;
        double km2;
        double knew;
        while (true) {
            //cout<<"cont2: "<<cont2<<endl;
            /*
             % Begin block 1
             %
             % Compute coefficients of formulas for this step. Avoid computing
             % those quantities not changed when step size is not changed.
             */
            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;

            //ns is the number of steps taken with size h, including the
            //current one. When k<ns, no coefficients change.

            if (h != hold) {
                ns = 0;
            }
            if (ns <= kold) {
                ns = ns + 1;
            }
            int nsp1 = ns + 1;

            if (k >= ns) {
                //Compute those components of alpha[*],beta[*],psi[*],sig[*]
                //which are changed
                beta(ns, 0) = 1.0;
                int realns = ns;
                alpha(ns, 0) = 1.0 / realns;
                double temp1 = h * realns;
                sig(nsp1, 0) = 1.0;
                if (k >= nsp1) {
                    for (int i = nsp1; i <= k; i++) {
                        int im1 = i - 1;
                        double temp2 = psi_(im1, 0);
                        psi_(im1, 0) = temp1;
                        beta(i, 0) = beta(im1, 0) * psi_(im1, 0) / temp2;
                        temp1 = temp2 + h;
                        alpha(i, 0) = h / temp1;
                        int reali = i;
                        sig(i + 1, 0) = reali * alpha(i, 0) * sig(i, 0);
                    }
                }
                psi_(k, 0) = temp1;

                //Compute coefficients g[*]; initialize v[*] and set w[*].
                if (ns > 1) {
                    //If order was raised, update diagonal part of v[*]
                    if (k > kold) {
                        double temp4 = k * kp1;
                        v(k, 0) = 1.0 / temp4;
                        int nsm2 = ns - 2;
                        for (int j = 0; j < nsm2; j++) {
                            int i = k - j - 2;
                            v(i + 1, 0) = v(i + 1, 0) - alpha(j + 2, 0) * v(i + 2, 0);
                        }
                    }

                    //Update V[*] and set W[*]
                    double limit1 = kp1 - ns;
                    double temp5 = alpha(ns, 0);
                    for (int iq = 0; iq < limit1; iq++) {
                        v(iq + 1, 0) = v(iq + 1, 0) - temp5 * v(iq + 2, 0);
                        w(iq + 1, 0) = v(iq + 1, 0);
                    }
                    g(nsp1, 0) = w(1, 0);

                } else {
                    for (int iq = 1; iq <= k; iq++) {
                        double temp3 = iq * (iq + 1);
                        v(iq, 0) = 1.0 / temp3;
                        w(iq, 0) = v(iq, 0);
                    }
                }

                //Compute the g[*] in the work vector w[*]
                int nsp2 = ns + 2;
                if (kp1 >= nsp2) {
                    for (int i = nsp2; i <= kp1; i++) {
                        double limit2 = kp2 - i;
                        double temp6 = alpha(i - 1, 0);
                        for (int iq = 0; iq < limit2; iq++) {
                            w(iq + 1, 0) = w(iq + 1, 0) - temp6 * w(iq + 2, 0);
                        }
                        //w.print();
                        g(i, 0) = w(1, 0);
                    }
                }
            } // if K>=NS

            /*
             * End block 1
             */

            /*
             * Begin block 2
            %
            % Predict a solution p[*], evaluate derivatives using predicted
            % solution, estimate local error at order k and errors at orders
            % k, k-1, k-2 as if constant step size were used.
             */

            // Change phi to phi star
            if (k >= nsp1) {
                for (int i = nsp1; i <= k; i++) {
                    double temp1 = beta(i, 0);
                    for (int l = 0; l < n_eqn; l++) {
                        phi(l, i ) = temp1 * phi(l, i );
                    }
                }
            }

            //Predict solution and differences
            for (int l = 0; l < n_eqn; l++) {
                phi(l, kp2) = phi(l, kp1);
                phi(l, kp1 ) = 0.0;
                p(l, 0) = 0.0;
            }
            for (int j = 1; j <= k; j++) {
                int i = kp1 - j;
                int ip1 = i + 1;
                double temp2 = g(i , 0);
                for (int l = 0; l < n_eqn; l++) {
                    p(l, 0) = p(l, 0) + temp2 * phi(l, i );
                    phi(l, i ) = phi(l, i ) + phi(l, ip1);
                }
            }
            if (nornd) {
                //y.print();
                //p.print();

                p = y + (p*h);
            } else {
                for (int l = 0; l < n_eqn; l++) {
                    double tau = h * p(l, 0) - phi(l, 15);
                    p(l, 0) = y(l, 0) + tau;
                    phi(l, 16) = (p(l, 0) - y(l, 0)) - tau;
                }
            }
            double xold = x;
            x = x + h;
            absh = fabs(h);
            yp = func(x, p);

            // Estimate errors at orders k, k-1, k-2
            double erkm2 = 0.0;
            erkm1 = 0.0;
            double erk = 0.0;

            for (int l = 0; l < n_eqn; l++) {
                double temp3 = 1.0 / wt(l, 0);
                //yp.print();
                //phi.print();
                //cout<<setprecision(40)<<yp(l,0)<<endl;
                //cout<<setprecision(40)<<phi(l, 1)<<endl;
                double temp4 = yp(l, 0) - phi(l, 1);
                if (km2 > 0.0) {
                    erkm2 = erkm2 + ((phi(l, km1) + temp4) * temp3) * ((phi(l, km1) + temp4) * temp3);
                }
                if (km2 >= 0.0) {
                    erkm1 = erkm1 + ((phi(l, k) + temp4) * temp3) * ((phi(l, k) + temp4) * temp3);
                }
                erk = erk + (temp4 * temp3) * (temp4 * temp3);
            }
            if (km2 > 0) {
                erkm2 = absh * sig(km1, 0) * gstr(0,km2) * sqrt(erkm2);
            }
            if (km2 >= 0) {
                erkm1 = absh * sig(k, 0) * gstr(0,km1) * sqrt(erkm1);
            }
            double temp5 = absh * sqrt(erk);
            //g.print();
            double err = temp5 * (g(k, 0) - g(kp1, 0));
            //sig.print();
            //gstr.print();
            erk = temp5 * sig(kp1, 0) * gstr(0,k);
            knew = k;

            //Test if order should be lowered
            if (km2 > 0) {
                if (max(erkm1, erkm2) <= erk) {
                    knew = km1;
                }
            }
            if (km2 == 0) {
                if (erkm1 <= 0.5 * erk) {
                    knew = km1;
                }
            }

            /*
             * End block 2
             */

            /*
             * If step is successful continue with block 4, otherwise repeat
             * blocks 1 and 2 after executing block 3
             */
            /*if(cont2>=19){
                err=err*0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001;
            }*/
            bool success = (err <= epsilon);
            cout<<"cont2: "<<cont2<<endl;
            cont2++;
            if (!success) {
                /*
                 * Begin block 3
                 */

                /*
                 * The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
                % 3rd consecutive failure, set order to 1. If step fails more
                % than 3 times, consider an optimal step size. Double error
                % tolerance and return if estimated step size is too small
                % for machine precision.
                 */

                // Restore x, phi[*,*] and psi[*]
                bool phase1 = false;
                x = xold;
                for (int i = 1; i <= k; i++) {
                    double temp1 = 1.0 / beta(i, 0);
                    int ip1 = i + 1;
                    for (int l = 0; l < n_eqn; i++) {
                        phi(l, i) = temp1 * (phi(l, i ) - phi(l, ip1));
                    }
                }

                if (k >= 2) {
                    for (int i = 2; i <= k; i++) {
                        psi_(i - 1, 0) = psi_(i, 0) - h;
                    }
                }

                // On third failure, set order to one.
                // Thereafter, use optimal step size
                ifail = ifail + 1;
                double temp2 = 0.5;
                if (ifail > 3) {
                    if (p5eps < 0.25 * erk) {
                        temp2 = sqrt(p5eps / erk);
                    }
                }
                if (ifail >= 3) {
                    knew = 1;
                }
                h = temp2 * h;
                k = knew;
                if (fabs(h) < fouru * fabs(x)) {
                    crash = true;
                    h = sign_(fouru * fabs(x), h);
                    epsilon = epsilon * 2.0;
                    return y;         //Exit
                }

                /*
                 * End block 3, return to start of block 1
                 */

            }  //end if(success)

            if (success) {
                break;
            }

        }

        /*
         * Begin block 4
         */

        /*
         * The step is successful. Correct the predicted solution, evaluate
         * the derivatives using the corrected solution and update the
         * differences. Determine best order and step size for next step.
         */

        kold = k;
        hold = h;

        //Correct and evaluate
        double temp1 = h * g(kp1, 0);
        if (nornd) {
            for (int l = 0; l < n_eqn; l++) {
                y(l, 0) = p(l, 0) + temp1 * (yp(l, 0) - phi(l, 1));
            }
        } else {
            rho = Matrix(1, 1);
            for (int l = 0; l < n_eqn; l++) {
                rho(0, 0) = ((yp(l, 0) - phi(l, 1)) * temp1) - phi(l, 16);
                y(l, 0) = p(l, 0) + rho(0, 0);
                phi(l, 15) = (y(l, 0) - p(l, 0)) - rho(0, 0);
            }
        }
        yp = func(x, y);
        //yp.print();
        //Update differences for next step
        for (int l = 0; l < n_eqn; l++) {
            phi(l, kp1) = yp(l, 0) - phi(l, 1);
            phi(l, kp2) = phi(l, kp1) - phi(l, kp2);
        }
        for (int i = 0; i < k; i++) {
            for (int l = 0; l < n_eqn; l++) {
                phi(l, i + 1) = phi(l, i + 1) + phi(l, kp1);
            }
        }
        //phi.print();

        /*
         * Estimate error at order k+1 unless
         * - in first phase when always raise order,
         * - already decided to lower order,
         * - step size not constant so estimate unreliable
         */
        double erkp1 = 0.0;
        bool phase1;
        double erk;
        if ((knew == km1) || (k == 12)) {
            phase1 = false;
        }
        if (phase1) {
            k = kp1;
            erk = erkp1;
        } else {
            if (knew == km1) {
                //lower order
                k = km1;
                erk = erkm1;
            } else {
                if (kp1 <= ns) {
                    for (int l = 0; l < n_eqn; l++) {
                        erkp1 = erkp1 + (phi(l, kp2) / wt(l, 0)) * (phi(l, kp2) / wt(l, 0));
                    }
                    erkp1 = absh * gstr(0,kp1) * sqrt(erkp1);
                    //Using estimated error at order k+1, determine
                    //appropriate order for next step
                    if (k > 1) {
                        if (erkm1 <= min(erk, erkp1)) {
                            //lower order
                            k = km1;
                            erk = erkm1;
                        } else {
                            if ((erkp1 < erk) && (k != 12)) {
                                //raise order
                                k = kp1;
                                erk = erkp1;
                            }
                        }
                    } else if (erkp1 < 0.5 * erk) {
                        /*
                         * raise order
                        % Here erkp1 < erk < max(erkm1,ermk2) else
                        % order would have been lowered in block 2.
                        % Thus order is to be raised
                         */
                        k = kp1;
                        erk = erkp1;
                    }
                }  // end if kp1<=ns
            }   //    end if knew!=km1
        }   //   end if !phase1

        //With new order determine appropriate step size for next step
        if (phase1 || (p5eps >= erk * two(0,k + 1))) {
            hnew = 2.0 * h;
        } else {
            if (p5eps < erk) {
                double temp2 = k + 1;
                double r = p5eps / pow(erk, 1.0 / temp2);
                hnew = absh * max(0.5, min(0.9, r));
                hnew = sign_(max(hnew, fouru * fabs(x)), h);
            } else {
                hnew = h;
            }
        }
        h = hnew;

        /*
         * End block 4
         */


        //Test for too small tolerances
        if (crash) {
            State_ = DE_STATE.DE_BADACC;
            relerr = epsilon * releps;       // Modify relative and absolute
            abserr = epsilon * abseps;       // accuracy requirements
            y = yy;                   // Copy last step
            t = x;
            told = t;
            bool OldPermit = true;
            return y;                       // Weak failure exit
        }

        nostep = nostep + 1;  // Count total number of steps

        //Count number of consecutive steps taken with the order of
        //the method being less or equal to four and test for stiffness

        kle4 = kle4 + 1;
        if (kold > 4) {
            kle4 = 0;
        }
        if (kle4 >= 50) {
            stiff = true;
        }
    }
    cont ++;
    return Matrix(0,0); //No se llegará nunca aquí luego podemos devolver lo que sea.
} //End step loop
/*
%if ( State_==DE_STATE.DE_INVPARAM )
%       error ('invalid parameters in DEInteg');
%       exit;
%   end
%   if ( State_==DE_STATE.DE_BADACC )
%       warning ('on','Accuracy requirement not achieved in DEInteg');
%   end
%   if ( State_==DE_STATE.DE_STIFF )
%       warning ('on','Stiff problem suspected in DEInteg');
%   end
%   if ( State_ >= DE_STATE.DE_DONE )
%       break;
%   end
%
% end
 */