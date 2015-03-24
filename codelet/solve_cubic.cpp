/*
http://tog.acm.org/resources/GraphicsGems/
EULA: The Graphics Gems code is copyright-protected. In other words,
you cannot claim the text of the code as your own and resell it.
Using the code is permitted in any program, product, or library,
non-commercial or commercial. Giving credit is not required,
though is a nice gesture. The code comes as-is, and if there are
any flaws or problems with any Gems code, nobody involved with Gems
- authors, editors, publishers, or webmasters - are to be held responsible.
 Basically, don't be a jerk, and remember that anything free comes 
 with no guarantee.
*/

/********************************************************
*                                                       *
* This function determines the roots of a cubic         *
* equation.                                             *
* It takes as parameters a pointer to the four          *
* coefficient of the cubic equation (the c[3] is the    *
* coefficient of x3 and so on) and a pointer to the     *
* three element array in which the roots are to be      *
* placed.                                               *
* It outputs the number of roots found                  *
*                                                       *
********************************************************/

int solveCubic(double c[4], double s[3])
{
int     i, num;
double  sub,
        A, B, C,
        sq_A, p, q,
        cb_p, D;

// normalize the equation:x ^ 3 + Ax ^ 2 + Bx  + C = 0
A = c[2] / c[3];
B = c[1] / c[3];
C = c[0] / c[3];

// substitute x = y - A / 3 to eliminate the quadric term: x^3 + px + q = 0

sq_A = A * A;
p = 1.0/3.0 * (-1.0/3.0 * sq_A + B);
q = 1.0/2.0 * (2.0/27.0 * A *sq_A - 1.0/3.0 * A * B + C);

// use Cardano's formula

cb_p = p * p * p;
D = q * q + cb_p;

if (D == 0)
    {
    if (q == 0)
        {
        // one triple solution
        s[0] = 0.0;
        num = 1;
        }
    else
        {
        // one single and one double solution
        double u = cbrt(-q);
        s[0] = 2.0 * u;
        s[1] = - u;
        num = 2;
        }
    }
else
    if (D < 0.0)
        {
        // casus irreductibilis: three real solutions
        double phi = 1.0/3.0 * acos(-q / sqrt(-cb_p));
        double t = 2.0 * sqrt(-p);
        s[0] = t * cos(phi);
        s[1] = -t * cos(phi + M_PI / 3.0);
        s[2] = -t * cos(phi - M_PI / 3.0);
        num = 3;
        }
    else
        {
        // one real solution
        double sqrt_D = sqrt(D);
        double u = cbrt(sqrt_D + fabs(q));
        if (q > 0.0)
            s[0] = - u + p / u ;
        else
            s[0] = u - p / u;
        num = 1;
        }

// resubstitute
sub = 1.0 / 3.0 * A;
for (i = 0; i < num; i++)
    s[i] -= sub;
return num;
}

