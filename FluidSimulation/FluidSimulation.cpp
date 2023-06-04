// FluidSimulation.cpp : Defines the entry point for the application.
//
// https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/GDC03.pdf
// https://mikeash.com/pyblog/fluid-simulation-for-dummies.html


#include "framework.h"
#include "FluidSimulation.h"
#include <windowsx.h>
#include <math.h>

#define MAX_LOADSTRING      100

#define CUBE_MATRIX_N       32
#define FLUID_BOX_SIZE      10
#define SIMULATION_DT       0.1f
#define WINDOW_WIDTH        (CUBE_MATRIX_N * FLUID_BOX_SIZE)
#define WINDOW_HEIGHT       (CUBE_MATRIX_N * FLUID_BOX_SIZE)

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name

// Forward declarations of functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)

struct FluidCube {
    int size;
    float dt;
    float diff;
    float visc;

    float* s;
    float* density;

    float* Vx;
    float* Vy;
    float* Vz;

    float* Vx0;
    float* Vy0;
    float* Vz0;
};

FluidCube* CubeMatrix;

FluidCube* FluidCubeCreate(int size, float diffusion, float viscosity, float dt)
{
    FluidCube* cube = (FluidCube *)malloc(sizeof(*cube));
    int N = size;

    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;

    cube->s = (float *)calloc(N * N * N, sizeof(float));
    cube->density = (float*)calloc(N * N * N, sizeof(float));

    cube->Vx = (float*)calloc(N * N * N, sizeof(float));
    cube->Vy = (float*)calloc(N * N * N, sizeof(float));
    cube->Vz = (float*)calloc(N * N * N, sizeof(float));

    cube->Vx0 = (float*)calloc(N * N * N, sizeof(float));
    cube->Vy0 = (float*)calloc(N * N * N, sizeof(float));
    cube->Vz0 = (float*)calloc(N * N * N, sizeof(float));

    return cube;
}

void FluidCubeFree(FluidCube* cube)
{
    free(cube->s);
    free(cube->density);

    free(cube->Vx);
    free(cube->Vy);
    free(cube->Vz);

    free(cube->Vx0);
    free(cube->Vy0);
    free(cube->Vz0);

    free(cube);
}


static void set_bnd(int b, float* x, int N)
{
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
            x[IX(i, j, N - 1)] = b == 3 ? -x[IX(i, j, N - 2)] : x[IX(i, j, N - 2)];
        }
    }
    for (int k = 1; k < N - 1; k++) {
        for (int i = 1; i < N - 1; i++) {
            x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
            x[IX(i, N - 1, k)] = b == 2 ? -x[IX(i, N - 2, k)] : x[IX(i, N - 2, k)];
        }
    }
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
            x[IX(N - 1, j, k)] = b == 1 ? -x[IX(N - 2, j, k)] : x[IX(N - 2, j, k)];
        }
    }

    x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)]
        + x[IX(0, 1, 0)]
        + x[IX(0, 0, 1)]);
    x[IX(0, N - 1, 0)] = 0.33f * (x[IX(1, N - 1, 0)]
        + x[IX(0, N - 2, 0)]
        + x[IX(0, N - 1, 1)]);
    x[IX(0, 0, N - 1)] = 0.33f * (x[IX(1, 0, N - 1)]
        + x[IX(0, 1, N - 1)]
        + x[IX(0, 0, N)]);
    x[IX(0, N - 1, N - 1)] = 0.33f * (x[IX(1, N - 1, N - 1)]
        + x[IX(0, N - 2, N - 1)]
        + x[IX(0, N - 1, N - 2)]);
    x[IX(N - 1, 0, 0)] = 0.33f * (x[IX(N - 2, 0, 0)]
        + x[IX(N - 1, 1, 0)]
        + x[IX(N - 1, 0, 1)]);
    x[IX(N - 1, N - 1, 0)] = 0.33f * (x[IX(N - 2, N - 1, 0)]
        + x[IX(N - 1, N - 2, 0)]
        + x[IX(N - 1, N - 1, 1)]);
    x[IX(N - 1, 0, N - 1)] = 0.33f * (x[IX(N - 2, 0, N - 1)]
        + x[IX(N - 1, 1, N - 1)]
        + x[IX(N - 1, 0, N - 2)]);
    x[IX(N - 1, N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1, N - 1)]
        + x[IX(N - 1, N - 2, N - 1)]
        + x[IX(N - 1, N - 1, N - 2)]);
}

static void lin_solve(int b, float* x, float* x0, float a, float c, int iter, int N)
{
    float cRecip = 1.0f / c;
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m)] =
                        (x0[IX(i, j, m)]
                            + a * (x[IX(i + 1, j, m)]
                                + x[IX(i - 1, j, m)]
                                + x[IX(i, j + 1, m)]
                                + x[IX(i, j - 1, m)]
                                + x[IX(i, j, m + 1)]
                                + x[IX(i, j, m - 1)]
                                )) * cRecip;
                }
            }
        }
        set_bnd(b, x, N);
    }
}

static void diffuse(int b, float* x, float* x0, float diff, float dt, int iter, int N)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void advect(int b, float* d, float* d0, float* velocX, float* velocY, float* velocZ, float dt, int N)
{
    float i0, i1, j0, j1, k0, k1;

    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    float dtz = dt * (N - 2);

    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y, z;

    float Nfloat = (float)N;
    float ifloat, jfloat, kfloat;
    int i, j, k;

    for (k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
        for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
            for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j, k)];
                tmp2 = dty * velocY[IX(i, j, k)];
                tmp3 = dtz * velocZ[IX(i, j, k)];
                x = ifloat - tmp1;
                y = jfloat - tmp2;
                z = kfloat - tmp3;

                if (x < 0.5f) x = 0.5f;
                if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
                i0 = floorf(x);
                i1 = i0 + 1.0f;
                if (y < 0.5f) y = 0.5f;
                if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
                j0 = floorf(y);
                j1 = j0 + 1.0f;
                if (z < 0.5f) z = 0.5f;
                if (z > Nfloat + 0.5f) z = Nfloat + 0.5f;
                k0 = floorf(z);
                k1 = k0 + 1.0f;

                s1 = x - i0;
                s0 = 1.0f - s1;
                t1 = y - j0;
                t0 = 1.0f - t1;
                u1 = z - k0;
                u0 = 1.0f - u1;

                int i0i = (int)i0;
                int i1i = (int)i1;
                int j0i = (int)j0;
                int j1i = (int)j1;
                int k0i = (int)k0;
                int k1i = (int)k1;

                if (i0i >= N) i0i = (N - 1);
                if (i1i >= N) i1i = (N - 1);
                if (j0i >= N) j0i = (N - 1);
                if (j1i >= N) j1i = (N - 1);
                if (k0i >= N) k0i = (N - 1);
                if (k1i >= N) k1i = (N - 1);

                d[IX(i, j, k)] =

                    s0 * (t0 * (u0 * d0[IX(i0i, j0i, k0i)]
                        + u1 * d0[IX(i0i, j0i, k1i)])
                        + (t1 * (u0 * d0[IX(i0i, j1i, k0i)]
                            + u1 * d0[IX(i0i, j1i, k1i)])))
                    + s1 * (t0 * (u0 * d0[IX(i1i, j0i, k0i)]
                        + u1 * d0[IX(i1i, j0i, k1i)])
                        + (t1 * (u0 * d0[IX(i1i, j1i, k0i)]
                            + u1 * d0[IX(i1i, j1i, k1i)])));
            }
        }
    }
    set_bnd(b, d, N);
}

static void project(float* velocX, float* velocY, float* velocZ, float* p, float* div, int iter, int N)
{
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j, k)] = -0.5f * (
                    velocX[IX(i + 1, j, k)]
                    - velocX[IX(i - 1, j, k)]
                    + velocY[IX(i, j + 1, k)]
                    - velocY[IX(i, j - 1, k)]
                    + velocZ[IX(i, j, k + 1)]
                    - velocZ[IX(i, j, k - 1)]
                    ) / N;
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
                velocX[IX(i, j, k)] -= 0.5f * (p[IX(i + 1, j, k)]
                    - p[IX(i - 1, j, k)]) * N;
                velocY[IX(i, j, k)] -= 0.5f * (p[IX(i, j + 1, k)]
                    - p[IX(i, j - 1, k)]) * N;
                velocZ[IX(i, j, k)] -= 0.5f * (p[IX(i, j, k + 1)]
                    - p[IX(i, j, k - 1)]) * N;
            }
        }
    }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
    set_bnd(3, velocZ, N);
}

void FluidCubeStep(FluidCube* cube)
{
    int N = cube->size;
    float visc = cube->visc;
    float diff = cube->diff;
    float dt = cube->dt;
    float* Vx = cube->Vx;
    float* Vy = cube->Vy;
    float* Vz = cube->Vz;
    float* Vx0 = cube->Vx0;
    float* Vy0 = cube->Vy0;
    float* Vz0 = cube->Vz0;
    float* s = cube->s;
    float* density = cube->density;

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

void FluidCubeAddDensity(FluidCube* cube, int x, int y, int z, float amount)
{
    int N = cube->size;
    cube->density[IX(x, y, z)] += amount;
}

void FluidCubeAddVelocity(FluidCube* cube, int x, int y, int z, float amountX, float amountY, float amountZ)
{
    int N = cube->size;
    int index = IX(x, y, z);

    cube->Vx[index] += amountX;
    cube->Vy[index] += amountY;
    cube->Vz[index] += amountZ;
}

void DrawBox(HDC hdc, int x, int y, float density)
{
    RECT rect = { 0 };
    UINT grayScale = (UINT)(255.0f * density);
    grayScale = (grayScale >= 255 ? 255 : grayScale);
    HBRUSH brush = CreateSolidBrush(RGB(grayScale, grayScale, grayScale));

    int px = (x * FLUID_BOX_SIZE);
    int py = (y * FLUID_BOX_SIZE);

    SetRect(&rect, px, py, (px + FLUID_BOX_SIZE), (py + FLUID_BOX_SIZE));
    FillRect(hdc, &rect, brush);
    DeleteObject(brush);
}

void SetClientSize(HWND hwnd, int clientWidth, int clientHeight)
{
    if (IsWindow(hwnd))
    {

        DWORD dwStyle = (DWORD)GetWindowLongPtr(hwnd, GWL_STYLE);
        DWORD dwExStyle = (DWORD)GetWindowLongPtr(hwnd, GWL_EXSTYLE);
        HMENU menu = GetMenu(hwnd);

        RECT rc = { 0, 0, clientWidth, clientHeight };

        AdjustWindowRectEx(&rc, dwStyle, menu ? TRUE : FALSE, dwExStyle);
            
        SetWindowPos(hwnd, NULL, 0, 0, rc.right - rc.left, rc.bottom - rc.top,
            SWP_NOZORDER | SWP_NOMOVE);
    }
}

VOID CALLBACK MyTimerProc(
    HWND hWnd,        // handle to window for timer messages 
    UINT message,     // WM_TIMER message 
    UINT idTimer,     // timer identifier 
    DWORD dwTime)     // current system time 
{
    FluidCubeStep(CubeMatrix);
    InvalidateRect(hWnd, NULL, FALSE);
}

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    CubeMatrix = FluidCubeCreate(CUBE_MATRIX_N, 0.0001f, 0.001f, SIMULATION_DT);

    // Initialize global strings
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_FLUIDSIMULATION, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Perform application initialization:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_FLUIDSIMULATION));

    MSG msg;

    // Main message loop:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}



//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_FLUIDSIMULATION));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)CreateSolidBrush(0x00000000);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_FLUIDSIMULATION);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // Store instance handle in our global variable

   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
       CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   SetClientSize(hWnd, WINDOW_WIDTH, WINDOW_HEIGHT);
   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   SetTimer(hWnd, 1, (UINT)(SIMULATION_DT * 1000), (TIMERPROC)MyTimerProc);

   return TRUE;
}

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE: Processes messages for the main window.
//
//  WM_COMMAND  - process the application menu
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    static int xPosPrev = 0;
    static int yPosPrev = 0;

    switch (message)
    {
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // Parse the menu selections:
            switch (wmId)
            {
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            int N = CubeMatrix->size;

            for (int i = 0, ix = 0, iy = 0; i < (N * N); i++, ix++)
            {
                if (ix >= N)
                {
                    ix = 0;
                    iy++;
                }

                DrawBox(hdc, ix, iy, CubeMatrix->density[IX(ix, iy, (N/2))]);
            }
            EndPaint(hWnd, &ps);
        }
        break;
    case WM_MOUSEMOVE:
        {
            if (wParam & MK_LBUTTON)
            {
                int xPos = GET_X_LPARAM(lParam);
                int yPos = GET_Y_LPARAM(lParam);
                FluidCubeAddDensity(CubeMatrix, (xPos / FLUID_BOX_SIZE), (yPos / FLUID_BOX_SIZE), (CUBE_MATRIX_N / 2), 1.0f);
                FluidCubeAddVelocity(CubeMatrix, (xPos / FLUID_BOX_SIZE), (yPos / FLUID_BOX_SIZE), (CUBE_MATRIX_N / 2), (float)(xPos - xPosPrev), (float)(yPos - yPosPrev), 0);
                xPosPrev = xPos;
                yPosPrev = yPos;
            }
        }
        break;
    case WM_DESTROY:
        KillTimer(hWnd, 1);
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
