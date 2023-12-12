#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <dir.h>
#include "libraries/pbPlots.c"
#include "libraries/supportLib.c"

#define PI 3.14159f

// Global Variables
char rules[2][2];
char rules_list[3] = {'N', 'Z', 'P'};
double mu[2][3];

double sgn(double x){
    if(x>0){
        return 1;
    }else if(x<0){
        return -1;
    }else{
        return 0;
    }
}

// Function to add char in a string
char add(char *s, const int len, char x, int pos){
    s[len] = 0;
    memmove(s+pos+1, s+pos, len-pos+1);
    s[pos]=x;
}

// Function to find the index from rules_list
int find_idx(char x){
    int idx;
    for(idx = 0; idx<5; idx++){
        if(rules_list[idx] == x){
            break;
        }
    }

    return idx;
}

// Function to calculate mu
double *calculate_mu(double error, double a, double b, int idx){
    static double mu_result[2];
    mu_result[0] = (b - error)/(b-a);
    mu_result[1] = (error - a)/(b-a);

    return mu_result;
}

// Find rules and calculate mu
void membership_function(double error, double *e, int idx){
    double e1 = e[0];
    double e2 = e[1];
    double e3 = e[2];
    double *mu_result;

    if(error <= e1){
        mu[idx][0] = 1;
        mu[idx][1] = 0;
        mu[idx][2] = 0;
        rules[idx][0] = 'N';
        rules[idx][1] = 'Z';
    }else if(error < e2){
        mu_result = calculate_mu(error, e1, e2, idx);
        mu[idx][0] = mu_result[0];
        mu[idx][1] = mu_result[1];
        mu[idx][2] = 0;
        rules[idx][0] = 'N';
        rules[idx][1] = 'Z';
    }else if(error == e2){
        mu[idx][0] = 0;
        mu[idx][1] = 1;
        mu[idx][2] = 0;
        rules[idx][0] = 'N';
        rules[idx][1] = 'Z';
    }else if(error < e3){
        mu_result = calculate_mu(error, e2, e3, idx);
        mu[idx][0] = 0;
        mu[idx][1] = mu_result[0];
        mu[idx][2] = mu_result[1];
        rules[idx][0] = 'Z';
        rules[idx][1] = 'P';
    }else if(error >= e3){
        mu[idx][0] = 0;
        mu[idx][1] = 0;
        mu[idx][2] = 1;
        rules[idx][0] = 'Z';
        rules[idx][1] = 'P';
    }
}

void fuzzy_control(double defuzz[], double err_range[], double delta_err_range[], int num_plant, int iteration, double lr, int epoch){
    // Output Variable
    double y[5000] = {0.0};
    double x2[5000] = {0.0};
    double x1[5000] = {0.0};
    
    // Fuzzy output variable
    double u[5000] = {0.0};
    // Desired Output Model
    double ym[5000] = {0.0};
    // Time variable
    double time[5000];
    // Previous errors
    double err_prev;
    double delta_err_prev;
    // Current errors
    double err = ym[0] - y[0];
    double delta_err;
    double double_delta_err;

    // Get all defuzzifier value
    double N = defuzz[0];
    double Z = defuzz[1];
    double P = defuzz[2];
    
    // Neural network variables
    double w[9][4]; // Wmn
    double R[3];    // Rm
    double z[9];    // Zm
    double net[4];  // Netn
    int R_idx[9];

    // Force of rules matrix
    double rules_matrix[3][3] = {{N, N, P},
                                 {N, Z, P},
                                 {N, P, P}};

    R[0] = N;
    R[1] = Z;
    R[2] = P;

    for(int m = 0; m<9; m++){
        for(int n = 0; n<4; n++){
            w[m][n] = 1;
        }
    }

    for(int m = 0; m<3; m++){
        for(int n = 0; n<3; n++){
            if(rules_matrix[m][n] == N){
                R_idx[(m*3)+n] = 0;
            }else if(rules_matrix[m][n] == Z){
                R_idx[(m*3)+n] = 1;
            }else if(rules_matrix[m][n] == P){
                R_idx[(m*3)+n] = 2;
            }
        }
    }

    // Set value for time and desired output
    for(int i = 0; i<iteration; i++){
        time[i] = (double)i;
    }

    for(int l = 0; l<epoch; l++){
        // Reset all variables
        int iter = 0;
        if(num_plant != 3){
            if(l == epoch - 1){
                iter = iteration;
            }else{
                iter = 500;
            }
        }else{
            iter = 2000;
        }
        for(int i = 0; i<iter; i++){
            // Output Variable
            y[i] = 0.0;
            x2[i] = 0.0;
            x1[i] = 0.0;
            
            // Fuzzy output variable
            u[i] = 0.0;
            // Desired Output Model
            ym[i] = 0.0;
        }
        err = ym[0] - y[0];

        for(int k = 0; k<iter-1; k++){
            // Error calculations
            err_prev = err;
            err = ym[k] - y[k];
            delta_err = err - err_prev;

            // printf("%f ", err);

            // Fuzzifier
            membership_function(err, err_range, 0);
            membership_function(delta_err, delta_err_range, 1);

            // Calculate weights
            for(int i = 0; i<3; i++){
                for(int j = 0; j<3; j++){
                    z[(i*3)+j] = mu[1][i]*mu[0][j];
                }
            }

            // Find fired index
            int fired_idx[4];
            for(int i = 0; i<2; i++){
                int i_idx = find_idx(rules[1][i]); 
                
                for(int j = 0; j<2; j++){
                    int j_idx = find_idx(rules[0][j]);

                    fired_idx[(i*2)+j] = (i_idx*3)+j_idx;
                }
            }

            u[k+1] = 0;

            for(int n = 0; n<4; n++){
                double sigma = 0;
                for(int i = 0; i<4; i++){
                    int m = fired_idx[i];
                
                    sigma += z[m]*w[m][n];
                }

                int m = fired_idx[n];

                net[n] = z[m]*w[m][n] / sigma;

                u[k+1] += net[n]*R[R_idx[m]];
            }

            // Fuzzy output as Δu(k)
            u[k+1] = u[k] + u[k+1];

            if(num_plant == 1){
                // Model Reference 1
                double r = 0.5*sin(0.007*k) + 2*cos(0.059*k);
                ym[k+1] = (0.99*ym[k] + 0.01*r);
            }else if(num_plant == 2){
                // Model Reference 2
                double uc = sin((2*PI*k)/250);
                double fuc = 0.6*sin(PI*uc) + 0.3*sin(3*PI*uc) + 0.1*sin(5*PI*uc);
                ym[k+1] = (0.3*ym[k] + 0.6*ym[k-1] + fuc);
            }else if(num_plant == 3){
                // Model Reference 3
                double r;
                if(k < 500){
                    r = sin((PI*k)/25);
                }else if(k < 1000){
                    r = 1;
                }else if(k < 1500){
                    r = -1;
                }else{
                    r = (0.3*sin((PI*k)/25)) + (0.4*sin((PI*k)/32)) + (0.3*sin((PI*k)/40));
                }
                ym[k+1] = ((0.6*ym[k]) + (0.2*ym[k-1]) + (0.1*r));
            }

            if(num_plant == 1){
                // 1st Plant
                x1[k+1] = (x1[k]+(0.01*x2[k])+(0.01*u[k+1]));
                x2[k+1] = ((0.1*x1[k])+(0.97*x2[k]));
                y[k+1] = (x1[k+1]);
            }else if(num_plant == 2){
                // 2nd Plant
                y[k+1] = ((0.3*y[k]) + (0.6*y[k-1]) + u[k+1]);
            }else if(num_plant == 3){
                // 3rd Plant
                y[k+1] = (0.35*((y[k]*y[k-1]*(y[k]+2.5))/(1+pow(y[k], 2)+pow(y[k-1], 2))) + (0.35*u[k+1]));
            }

            if(l != epoch-1){

                double d;
                d = (ym[k+1] - y[k+1]) * sgn((y[k+1] - y[k])/(u[k] - u[k-1]));

                for(int n = 0; n<4; n++){
                    int m = fired_idx[n];

                    R[R_idx[m]] -= lr*(d*net[n]);

                    if(R_idx[m] == 0)
                        R[R_idx[m]] = fmin(R[R_idx[m]], -0.1);
                    if(R_idx[m] == 2)
                        R[R_idx[m]] = fmax(R[R_idx[m]], 0.1);

                    if (R[R_idx[m]] > 20)
                        R[R_idx[m]] = 20;
                    else if(R[R_idx[m]] < -20)
                        R[R_idx[m]] = -20;
                }

                for(int n = 0; n<4; n++){
                    double sigma = 0;
                    for(int i = 0; i<4; i++){
                        int m = fired_idx[i];

                        sigma += z[m]*w[m][n];
                    }

                    for(int m = 0; m<9; m++)
                        w[m][n] -= lr*(d * R[R_idx[m]] * ((((z[m] * sigma * w[m][n]) - (pow(z[m], 2) * w[m][n]))) / pow(sigma, 2)));
                }
            }
        }
    }

    printf("\n      R[1]     R[2]     R[3]");
    printf("\nR: ");
    for(int i=0; i<3; i++)
        printf("%f ", R[i]);

    printf("\n\n          w[m][1]   w[m][2]   w[m][3]   w[m][4]");
    for(int i=0; i<9; i++){
        printf("\nw[%d][n]: ", i+1);
        for(int j=0; j<4; j++)
            printf("%.7f ", w[i][j]);
    }
    printf("\n");

    // Make directory for saving plot results
    mkdir("outputs");

    char path[100] = "outputs/plant.png";
    char plant[2];
    char plant_str[10] = "PLANT ";
    wchar_t plant_str2[10];

    sprintf(plant, "%d", num_plant);

    add(path, 100, *plant, 13);
    add(plant_str, 10, *plant, 6);

    mbstowcs(plant_str2, plant_str, 10);

    // Plot output
    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
	series->xs = time;
	series->xsLength = iteration;
	series->ys = y;
	series->ysLength = iteration;
	series->linearInterpolation = true;
	series->lineType = L"dotted";
	series->lineTypeLength = wcslen(series->lineType);
	series->lineThickness = 2;
	series->color = CreateRGBColor(1, 0, 0);

    // Plot model reference
    ScatterPlotSeries *series2 = GetDefaultScatterPlotSeriesSettings();
	series2->xs = time;
	series2->xsLength = iteration;
	series2->ys = ym;
	series2->ysLength = iteration;
	series2->linearInterpolation = true;
	series2->lineType = L"solid";
	series2->lineTypeLength = wcslen(series2->lineType);
	series2->lineThickness = 2;
	series2->color = CreateRGBColor(0, 0, 1);

    // Plot settings
	ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
	settings->width = 1000;
	settings->height = 400;
	settings->autoBoundaries = true;
	settings->autoPadding = true;
    settings->title = plant_str2;
    settings->titleLength = wcslen(settings->title);
    settings->xLabel = L"Y(k)";
    settings->xLabelLength = wcslen(settings->xLabel);
    settings->yLabel = L"Steps";
    settings->yLabelLength = wcslen(settings->yLabel);
	ScatterPlotSeries *s [] = {series2, series};
	settings->scatterPlotSeries = s;
	settings->scatterPlotSeriesLength = 2;

	RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
	DrawScatterPlotFromSettings(canvasReference, settings);

	size_t length;
	double *pngdata = ConvertToPNG(&length, canvasReference->image);
	WriteToFile(pngdata, length, path);
	DeleteImage(canvasReference->image);

    printf("\n[INFO] Plot Saved To: '%s'\n", path);
}

int main(){
    // Defuzzifier parameter
    double uk[3] = {-10, 0, 10};
    // Error fuzzifier parameter
    double ek[3] = {-1, 0, 1};
    // Delta error fuzzifier parameter
    double delta_ek[3] = {-1, 0, 1};
    // Learning rate
    double lr;

    // Calculate and Plot 1st Plant
    lr = 1e-4;
    fuzzy_control(uk, ek, delta_ek, 1, 1000, lr, 5);

    // Calculate and Plot 2nd Plant
    lr = 1e-4;
    fuzzy_control(uk, ek, delta_ek, 2, 1000, lr, 500);

    // Calculate and Plot 3rd Plant
    lr = 5e-5;
    fuzzy_control(uk, ek, delta_ek, 3, 2000, lr, 430);

    printf("\n");

    return 0;
}