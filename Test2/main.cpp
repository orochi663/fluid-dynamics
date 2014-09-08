// simple test of a source code found on the internet

#include <math.h>
#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>


#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)


struct FluidCube {
    int size;
    float dt;
    float diff;
    float visc;
    
    float *s;
    float *density;
    
    float *Vx;
    float *Vy;
    float *Vz;

    float *Vx0;
    float *Vy0;
    float *Vz0;
};



FluidCube *FluidCubeCreate(int size, float diffusion, float viscosity, float dt)
{
    FluidCube *cube = new FluidCube();
    int N = size;
    
    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;
    
    cube->s = new float[N * N * N];
    cube->density = new float[N * N * N];
    
    cube->Vx = new float[N * N * N];
    cube->Vy = new float[N * N * N];
    cube->Vz = new float[N * N * N];
    
    cube->Vx0 = new float[N * N * N];
    cube->Vy0 = new float[N * N * N];
    cube->Vz0 = new float[N * N * N];
    
    return cube;
}

void FluidCubeFree(FluidCube *cube)
{
    delete[] cube->s;
    delete[] cube->density;
    
    delete[] cube->Vx;
    delete[] cube->Vy;
    delete[] cube->Vz;
    
    delete[] cube->Vx0;
    delete[] cube->Vy0;
    delete[] cube->Vz0;
    
    delete cube;
}


static void set_bnd(int b, float *x, int N)
{
    for(int j = 1; j < N - 1; j++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, j, 0  )] = b == 3 ? -x[IX(i, j, 1  )] : x[IX(i, j, 1  )];
            x[IX(i, j, N-1)] = b == 3 ? -x[IX(i, j, N-2)] : x[IX(i, j, N-2)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, 0  , k)] = b == 2 ? -x[IX(i, 1  , k)] : x[IX(i, 1  , k)];
            x[IX(i, N-1, k)] = b == 2 ? -x[IX(i, N-2, k)] : x[IX(i, N-2, k)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int j = 1; j < N - 1; j++) {
            x[IX(0  , j, k)] = b == 1 ? -x[IX(1  , j, k)] : x[IX(1  , j, k)];
            x[IX(N-1, j, k)] = b == 1 ? -x[IX(N-2, j, k)] : x[IX(N-2, j, k)];
        }
    }
    
    x[IX(0, 0, 0)]       = 0.33f * (x[IX(1, 0, 0)]
                                  + x[IX(0, 1, 0)]
                                  + x[IX(0, 0, 1)]);
    x[IX(0, N-1, 0)]     = 0.33f * (x[IX(1, N-1, 0)]
                                  + x[IX(0, N-2, 0)]
                                  + x[IX(0, N-1, 1)]);
    x[IX(0, 0, N-1)]     = 0.33f * (x[IX(1, 0, N-1)]
                                  + x[IX(0, 1, N-1)]
                                  + x[IX(0, 0, N)]);
    x[IX(0, N-1, N-1)]   = 0.33f * (x[IX(1, N-1, N-1)]
                                  + x[IX(0, N-2, N-1)]
                                  + x[IX(0, N-1, N-2)]);
    x[IX(N-1, 0, 0)]     = 0.33f * (x[IX(N-2, 0, 0)]
                                  + x[IX(N-1, 1, 0)]
                                  + x[IX(N-1, 0, 1)]);
    x[IX(N-1, N-1, 0)]   = 0.33f * (x[IX(N-2, N-1, 0)]
                                  + x[IX(N-1, N-2, 0)]
                                  + x[IX(N-1, N-1, 1)]);
    x[IX(N-1, 0, N-1)]   = 0.33f * (x[IX(N-2, 0, N-1)]
                                  + x[IX(N-1, 1, N-1)]
                                  + x[IX(N-1, 0, N-2)]);
    x[IX(N-1, N-1, N-1)] = 0.33f * (x[IX(N-2, N-1, N-1)]
                                  + x[IX(N-1, N-2, N-1)]
                                  + x[IX(N-1, N-1, N-2)]);
}

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m)] =
                        (x0[IX(i, j, m)]
                            + a*(    x[IX(i+1, j  , m  )]
                                    +x[IX(i-1, j  , m  )]
                                    +x[IX(i  , j+1, m  )]
                                    +x[IX(i  , j-1, m  )]
                                    +x[IX(i  , j  , m+1)]
                                    +x[IX(i  , j  , m-1)]
                           )) * cRecip;
                }
            }
        }
        set_bnd(b, x, N);
    }
}

static void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float *velocZ, float dt, int N)
{
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    float dtz = dt * (N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y, z;
    
    float Nfloat = N;
    float ifloat, jfloat, kfloat;
    int i, j, k;
    
    for(k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j, k)];
                tmp2 = dty * velocY[IX(i, j, k)];
                tmp3 = dtz * velocZ[IX(i, j, k)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                z    = kfloat - tmp3;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floorf(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floorf(y);
                j1 = j0 + 1.0f; 
                if(z < 0.5f) z = 0.5f;
                if(z > Nfloat + 0.5f) z = Nfloat + 0.5f;
                k0 = floorf(z);
                k1 = k0 + 1.0f;
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                u1 = z - k0;
                u0 = 1.0f - u1;
                
                int i0i = i0;
                int i1i = i1;
                int j0i = j0;
                int j1i = j1;
                int k0i = k0;
                int k1i = k1;
                
                d[IX(i, j, k)] = 
                
                    s0 * ( t0 * (u0 * d0[IX(i0i, j0i, k0i)]
                                +u1 * d0[IX(i0i, j0i, k1i)])
                        +( t1 * (u0 * d0[IX(i0i, j1i, k0i)]
                                +u1 * d0[IX(i0i, j1i, k1i)])))
                   +s1 * ( t0 * (u0 * d0[IX(i1i, j0i, k0i)]
                                +u1 * d0[IX(i1i, j0i, k1i)])
                        +( t1 * (u0 * d0[IX(i1i, j1i, k0i)]
                                +u1 * d0[IX(i1i, j1i, k1i)])));
            }
        }
    }
    set_bnd(b, d, N);
}

static void project(float *velocX, float *velocY, float *velocZ, float *p, float *div, int iter, int N)
{
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j, k)] = -0.5f*(
                         velocX[IX(i+1, j  , k  )]
                        -velocX[IX(i-1, j  , k  )]
                        +velocY[IX(i  , j+1, k  )]
                        -velocY[IX(i  , j-1, k  )]
                        +velocZ[IX(i  , j  , k+1)]
                        -velocZ[IX(i  , j  , k-1)]
                    )/N;
                p[IX(i, j, k)] = 0;
            }
        }
    }
    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, iter, N);
    
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j, k)] -= 0.5f * (  p[IX(i+1, j, k)]
                                                -p[IX(i-1, j, k)]) * N;
                velocY[IX(i, j, k)] -= 0.5f * (  p[IX(i, j+1, k)]
                                                -p[IX(i, j-1, k)]) * N;
                velocZ[IX(i, j, k)] -= 0.5f * (  p[IX(i, j, k+1)]
                                                -p[IX(i, j, k-1)]) * N;
            }
        }
    }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
    set_bnd(3, velocZ, N);
}

void FluidCubeStep(FluidCube *cube)
{
    int N          = cube->size;
    float visc     = cube->visc;
    float diff     = cube->diff;
    float dt       = cube->dt;
    float *Vx      = cube->Vx;
    float *Vy      = cube->Vy;
    float *Vz      = cube->Vz;
    float *Vx0     = cube->Vx0;
    float *Vy0     = cube->Vy0;
    float *Vz0     = cube->Vz0;
    float *s       = cube->s;
    float *density = cube->density;
    
    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    diffuse(2, Vy0, Vy, visc, dt, 4, N);
    diffuse(3, Vz0, Vz, visc, dt, 4, N);
    
    project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
    
    advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
    advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
    
    project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
    
    diffuse(0, s, density, diff, dt, 4, N);
    advect(0, density, s, Vx, Vy, Vz, dt, N);
}

void FluidCubeAddDensity(FluidCube *cube, int x, int y, int z, float amount)
{
    int N = cube->size;
    cube->density[IX(x, y, z)] += amount;
}

void FluidCubeAddVelocity(FluidCube *cube, int x, int y, int z, float amountX, float amountY, float amountZ)
{
    int N = cube->size;
    int index = IX(x, y, z);
    
    cube->Vx[index] += amountX;
    cube->Vy[index] += amountY;
    cube->Vz[index] += amountZ;
}

cv::Mat FluidCubeDrawVelocity(FluidCube *cube, int z)
{
    int N = cube->size;
    cv::Mat lImage(cube->size, cube->size, CV_32FC3);
    
    for(int x = 0; x < cube->size; ++x)
    {
        for(int y = 0; y < cube->size; ++y)
        {
            cv::Vec3f lValue(
                fabs(cube->Vx[IX(x,y,z)]),
                fabs(cube->Vy[IX(x,y,z)]),
                fabs(cube->Vz[IX(x,y,z)]));
            lImage.at<cv::Vec3f>(y, x) = 100.0f * lValue;
        }    
    }
    
    return lImage;
}

cv::Mat FluidCubeDrawDensity(FluidCube *cube, int z)
{
    int N = cube->size;
    cv::Mat lImage(cube->size, cube->size, CV_32FC3);
    
    for(int x = 0; x < cube->size; ++x)
    {
        for(int y = 0; y < cube->size; ++y)
        {
            cv::Vec3f lValue(fabs(cube->density[IX(x,y,z)]), 0.0f, 0.0f);
            lImage.at<cv::Vec3f>(y, x) = 100.0f * lValue;
        }    
    }
    
    return lImage;
}

int main(int argc, char** argv)
{
    cv::Mat lPotentials = cv::imread("gradient_cercle.png");
    std::vector<cv::Mat> lChannels;
    cv::split(lPotentials, lChannels);
    
    cv::Mat lForceX, lForceY;
    float lKernel[3] = {-1.0f, 0.0f, +1.0f};
    cv::filter2D(lChannels[0], lForceX, CV_32F, cv::Mat(cv::Size(1,3), CV_32FC1, lKernel));
    cv::filter2D(lChannels[0], lForceY, CV_32F, cv::Mat(cv::Size(3,1), CV_32FC1, lKernel));

    const int size = 64;
    FluidCube* lCube = FluidCubeCreate(size, 0, 0, 0.01f);
    
    char filename[255];
    
    for(int i = 0; i < 5000; ++i)
    {
        std::cout << "Step: " << i << std::endl;
        
        float angle = ((float)i / 200.0f) * (2.0f * 3.14159f);
        
        FluidCubeAddDensity(lCube, size/2, size/2, size/2, 50.0f);
        
        //FluidCubeAddVelocity(lCube, 2*size/3, size/2, size/2, 50.0f * cosf(angle), 50.0f * sinf(angle), 0.0f);
        //FluidCubeAddVelocity(lCube, 1*size/3, size/2, size/2, 50.0f * cosf(-angle+3.14159f*0.5f), 50.0f * sinf(-angle+3.14159f*0.5f), 0.0f);
        
        float clCoeff = 0.005f;
        
        for(int x = 0; x < size; ++x)
            for(int y = 0; y < size; ++y)
            {
                FluidCubeAddVelocity(
                    lCube, 
                    x,
                    y, 
                    size/2, 
                    lForceX.at<float>(y,x) * clCoeff,
                    lForceY.at<float>(y,x) * clCoeff,
                    0.0f);
            }
        
        FluidCubeStep(lCube);
        
        cv::Mat lVelocity = FluidCubeDrawVelocity(lCube, size/2);
        sprintf(filename, "velocity/frame_%05d.png", i);
        cv::imwrite(filename, lVelocity);
    }
    
    FluidCubeFree(lCube);

    return 0;
}































