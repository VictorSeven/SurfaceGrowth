                                       ///================================================================//
                                        //*                           Victor Buendia's                    *//
                                        //*                            SURFACE GROWTH                     *//
                                        //*                           Simulation project                  *//
                                        //*                   for Critical Cooperative Phenomena          *//
                                        ///================================================================//

//Include all C++ standards
#include<iostream>
#include<cstdlib>
#include<random>
#include<fstream>
#include<vector>

//Defines for array sizes
#define XMAX 1000
#define YMAX 101000
#define MAXSURF 10000


using namespace std;

//Functions which write where every block is
void uniform_deposit(int its, vector< vector<uint8_t> >& tetris, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i);
void diffusion_deposit(int its, vector< vector<uint8_t> >& tetris, int h[XMAX], int surfaces,  vector<double>& correlation, mt19937& gen, uniform_int_distribution<int>& ran_i);
void ballistic_deposit(int its, vector< vector<uint8_t> >& tetris, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i);
//This ones do the same, but only compute the correlation
void uniform_deposit(int its, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i);
void diffusion_deposit(int its, int h[XMAX], int surfaces,  vector<double>& correlation, mt19937& gen, uniform_int_distribution<int>& ran_i);
void ballistic_deposit(int its, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i);
//KPZ simulation
void kpz_probability(int h[XMAX], int surfaces, vector<int>& time,  vector<double>& correlation, mt19937& gen, uniform_real_distribution<double>& ran_u);
void kpz_milshtein(double h[XMAX], int surfaces,  vector<double>& correlation, mt19937& gen, normal_distribution<double>& ran_g);

//Overloads for computing correlation when h is in column for or continuous using KPZ
double corr_perp(int h[XMAX]);
double corr_perp(double h[XMAX]);


//Main program
int main(void)
{
    unsigned int i,j; //Counters

    unsigned int n; //Number of iterations
    unsigned int surfaces; //Number of desired surfaces
    unsigned int aver; //Number of averages

    int h[XMAX]; //Height of every colurmn
    //double h[XMAX];

    //vector<vector<uint8_t> > tetris(XMAX, vector<uint8_t>(YMAX)); //Different from 0 if it contains a particle. 8 bits to increase the avaiable space and allow colors

    ofstream output; //To write things

    //Random number generator and distribution
    mt19937 gen(958431198);
    uniform_int_distribution<int> ran_i(0, XMAX-1);
    uniform_real_distribution<double> ran_u(0.0, 1.0);
    normal_distribution<double> ran_g(0.0, 1.0);

    //Number of its and surfaces
    n = 1e8;
    surfaces = YMAX-1000;
    aver = 50;

    vector<double> correlation(2*surfaces, 0.0);
    vector<int> time(surfaces, 0.0);

    //Run the process
    cout << "start loop" << endl;
    for (i=0; i < aver; i++)
    {
        //diffusion_deposit(n, h, surfaces, correlation, gen, ran_i);
        //ballistic_deposit(n, h, surfaces, correlation, gen, ran_i);
        kpz_probability(h, surfaces, time, correlation, gen, ran_u);
        //kpz_milshtein(h, surfaces, correlation, gen, ran_g);
        cout << i << endl;
    }



    cout << "computed" << endl;
    //Open file and write a plot of the result
    /*output.open("viewsurfaces.txt");
    for (i=0; i < XMAX; i++)
    {
        for (j=0; j < YMAX; j++)
        {
            output << int(tetris[i][j]) << " ";
        }
        output << endl;
    }
    output.close();*/


    output.open("corrkpz.txt");
    for (i=0; i < surfaces; i++)
    {
        correlation[2*i] /= (1.0*aver);
        correlation[2*i+1] /= (1.0*aver);
        output << i << " " << correlation[2*i] << " " << correlation[2*i+1] - correlation[2*i]*correlation[2*i] << endl;
        if (correlation[2*i] == 0 and i > 10)
        {
            break;
        }
    }
    output.close();


    return 0;
}

///================================================================================================================
///Compute correlations and so on and additionally write where every block is
///================================================================================================================
//Blocks falling randomly
void uniform_deposit(int its, vector< vector<uint8_t> >& tetris, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i)
{
    int i,j;
    int index;

    //Init the height and block matrices
    for (i=0; i < XMAX; i++)
    {
        h[i] = 0;
        for (j = 0;  j < YMAX; j++)
        {
            tetris[i][j] = 0;
        }
    }


    its /= surfaces; //We are going to do its * surfaces iterations


    //For every surface I want,
    for (j = 0; j < surfaces; j++)
    {
        //Make the process
        for (i=0; i < its; i++)
        {
            index = ran_i(gen); //Get a column
            h[index] += 1; //Increase height of selected column
            tetris[index][h[index]] = j + 1; //Include this into the block matrix
        }

        correlation[j] = corr_perp(h);
    }
    return;
}

//Adds a diffusion to the model
void diffusion_deposit(int its, vector< vector<uint8_t> > &tetris, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i)
{
    int i,j;
    int index;

    int index_back, index_forw; //Indices for neighbours

    //Init the height and block matrices
    for (i=0; i < XMAX; i++)
    {
        h[i] = 0;
        for (j = 0;  j < YMAX; j++)
        {
            tetris[i][j] = 0;
        }
    }

    its /= surfaces; //We are going to do its * surfaces iterations

    for (j = 0; j < surfaces; j++)
    {
        for (i=0; i < its; i++)
        {
            index = ran_i(gen); //Get a column
            //Get nearest neighbours
            index_back = index == 0 ? XMAX-1 : index - 1;
            index_forw = index == XMAX-1 ? 0 : index + 1;

            //If one of the neighbours has is lower than our column,
            //then make the selected column the neighbour
            if (h[index_back] < h[index])   index = index_back;
            else if (h[index_forw] < h[index])  index = index_forw;

            //Once we know where we are going to fall, update height and blocks
            h[index] += 1;
            tetris[index][h[index]] = j % 2 + 1;
        }

        correlation[j] = corr_perp(h);

    }

    return;
}

//Ballistic deposition
void ballistic_deposit(int its, vector< vector<uint8_t> > &tetris, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i)
{
    int i,j;
    int index;

    int index_back, index_forw;

    //Init the height and block matrices
    for (i=0; i < XMAX; i++)
    {
        h[i] = 0;
        for (j = 0;  j < YMAX; j++)
        {
            tetris[i][j] = 0;
        }
    }

    its /= surfaces; //We are going to do its * surfaces iterations

    for (j = 0; j < surfaces; j++)
    {
        for (i=0; i < its; i++)
        {
            //Get the index and nearest neighbours
            index = ran_i(gen);
            index_back = index == 0 ? XMAX-1 : index - 1;
            index_forw = index == XMAX-1 ? 0 : index + 1;

            //If one of the nearest columns is higher than our, then the block
            //gets sticky and the height of our column becomes this
            if (h[index_back] > h[index])   h[index] = h[index_back];
            else if (h[index_forw] > h[index])  h[index] = h[index_forw];
            else h[index] += 1;

            //Also, update block matrix
            tetris[index][h[index]] = j + 1;
        }

        correlation[j] = corr_perp(h);
    }
    return;
}

///================================================================================================================
///Same, ot writing variable tetris
///================================================================================================================


void uniform_deposit(int its, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i)
{
    int i,j;
    int index;
    double corr;

    //Init the height and block matrices
    for (i=0; i < XMAX; i++)
    {
        h[i] = 0;
    }


    its /= surfaces; //We are going to do its * surfaces iterations


    //For every surface I want,
    for (j = 0; j < surfaces; j++)
    {
        //Make the process
        for (i=0; i < its; i++)
        {
            index = ran_i(gen); //Get a column
            h[index] += 1; //Increase height of selected column
        }

        corr = corr_perp(h);
        correlation[2*j] += corr;
        correlation[2*j+1] += corr * corr;
    }
    return;
}

//Adds a diffusion to the model
void diffusion_deposit(int its, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i)
{
    int i,j;
    int index;
    double corr;

    int index_back, index_forw; //Indices for neighbours

    //Init the height and block matrices
    for (i=0; i < XMAX; i++)
    {
        h[i] = 0;
    }

    its /= surfaces; //We are going to do its * surfaces iterations

    for (j = 0; j < surfaces; j++)
    {
        for (i=0; i < its; i++)
        {
            index = ran_i(gen); //Get a column
            //Get nearest neighbours
            index_back = index == 0 ? XMAX-1 : index - 1;
            index_forw = index == XMAX-1 ? 0 : index + 1;

            //If one of the neighbours has is lower than our column,
            //then make the selected column the neighbour
            if (h[index_back] < h[index])   index = index_back;
            else if (h[index_forw] < h[index])  index = index_forw;

            //Once we know where we are going to fall, update height and blocks
            h[index] += 1;
        }

        corr = corr_perp(h);
        correlation[2*j] += corr;
        correlation[2*j+1] += corr * corr;

    }

    return;
}

//Ballistic deposition
void ballistic_deposit(int its, int h[XMAX], int surfaces, vector<double>& correlation,  mt19937& gen, uniform_int_distribution<int>& ran_i)
{
    int i,j;
    int index;
    double corr;

    int index_back, index_forw;

    //Init the height and block matrices
    for (i=0; i < XMAX; i++)
    {
        h[i] = 0;
    }

    its /= surfaces; //We are going to do its * surfaces iterations

    for (j = 0; j < surfaces; j++)
    {
        for (i=0; i < its; i++)
        {
            //Get the index and nearest neighbours
            index = ran_i(gen);
            index_back = index == 0 ? XMAX-1 : index - 1;
            index_forw = index == XMAX-1 ? 0 : index + 1;

            //If one of the nearest columns is higher than our, then the block
            //gets sticky and the height of our column becomes this
            if (h[index_back] > h[index])   h[index] = h[index_back];
            else if (h[index_forw] > h[index])  h[index] = h[index_forw];
            else h[index] += 1;

        }

        corr = corr_perp(h);
        correlation[2*j] += corr;
        correlation[2*j+1] += corr * corr;
    }
    return;
}


//Computes the correlation, for Milstein
double corr_perp(double h[XMAX])
{
    int i;

    double hmean = 0.0;
    double aux;
    double corr = 0.0;

    for (i=0; i < XMAX; i++)
    {
        hmean += h[i];
    }

    hmean /= (1.0 * XMAX);

    for (i=0; i < XMAX; i++)
    {
        aux = (h[i] - hmean);
        corr += aux * aux;
    }

    corr = sqrt(corr/(1.0*XMAX));

    return corr;
}

//Computes the correlation
double corr_perp(int h[XMAX])
{
    int i;

    double hmean = 0.0;
    double hmean2 = 0.0;
    double aux;
    double corr = 0.0;


    for (i=0; i < XMAX; i++)
    {
        hmean += h[i];
        hmean2 += h[i] * h[i];
    }

    hmean /= (1.0 * XMAX);
    hmean2 /= (1.0 * XMAX);

    /*for (i=0; i < XMAX; i++)
    {
        aux = (h[i] - hmean);
        corr += aux * aux;
    } */

    //corr = sqrt(corr/(1.0*XMAX));
    corr = sqrt(hmean2 - hmean * hmean);
    return corr;
}


///================================================================================================================
///Do the same for KPZ using the algorithm by Mattos et al or Milstein
///================================================================================================================

void kpz_probability(int h[XMAX], int surfaces, vector<int>& time,  vector<double>& correlation, mt19937& gen, uniform_real_distribution<double>& ran_u)
{
    int i,j,k,m,n,t;
    int index;
    double corr;
    double counter;



    double rho = 0.5;
    double kappa = 0.1;
    double kernel;

    int index_back, index_forw; //Indices for neighbours

    //Init the height and block matrices
    for (i=0; i < XMAX; i++)
    {
        h[i] = 0;
    }


    //For the desired amount of surfaces
    for (j = 0; j < surfaces; j++)
    {
        //For every block,
        for (i=0; i < XMAX; i++)
        {
            //Get nearest neighbours
            index_back = i == 0 ? XMAX-1 : i - 1;
            index_forw = i == XMAX-1 ? 0 : i + 1;


            //Compute our kernel
            kernel = (h[index_forw] - h[index_back]) * (h[index_forw] - h[index_back]) / 2.0 + h[index_back] + h[index_forw] - 2.0 * h[i];

            //With probability p, put a block here
            if (ran_u(gen) <= rho * exp(kappa * kernel)) h[i] += 1;


        }

        corr = corr_perp(h);
        correlation[2*j] += corr;
        correlation[2*j+1] += corr * corr;

    }

    return;

}


void kpz_milshtein(double h[XMAX], int surfaces,  vector<double>& correlation, mt19937& gen, normal_distribution<double>& ran_g)
{
    int i,j,k;
    int index;
    double corr;

    double coef = 0.01 / 1.0; //h / deltaX^2

    int index_back, index_forw; //Indices for neighbours

    //Init the height and block matrices
    for (i=0; i < XMAX; i++)
    {
        h[i] = 0;
    }


    //For the desired amount of surfaces
    for (j = 0; j < surfaces; j++)
    {

            //For every block,
            for (i=0; i < XMAX; i++)
            {
                //Get nearest neighbours
                index_back = i == 0 ? XMAX-1 : i - 1;
                index_forw = i == XMAX-1 ? 0 : i + 1;
                h[i] += coef * (h[index_back] + h[index_forw] - 2.0 * h[i]) + coef * (h[index_forw] - h[index_back]) * (h[index_forw] - h[index_back]) / 4.0 + sqrt(coef) * ran_g(gen);
            }



        corr = corr_perp(h);
        correlation[2*j] += corr;
        correlation[2*j+1] += corr * corr;

    }
}
