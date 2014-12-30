// 2.cpp: определяет точку входа для консольного приложения.
//

 
#include "stdafx.h"
 
#include <fstream>
#include <iostream>
#include <iomanip>
 
#include <conio.h>
#include <Windows.h>
 
 
using namespace std;
 
inline double PI()
{
        __asm fldpi
}
 
int convert_x(const long double & x)
{
        return x * 640 / 0.001;
}
 
int convert_y(const long double & y)
{
        return y * 320 / 0.01;
}
 
/**long double KoefD(
        const double & y,
        const double & zz)
{
        long double x;
        x = (3 * kb * y) / (2 * pow(zz, 2), pow(e, 2));
        x *= pow(kb * Te / PI() * ne * pow(e, 2), 0.5);
        return x;
}
 
const long double
 
De = KoefD(Te, z),
Di = KoefD(Ti, z),
Dh = KoefD(Th, z);**/
 
/* 22 */
 
/* long double KoefK(
        const double & a,
        const double & b,
        const double & c,
        const double & d)
{
        long double kk;
        kk = pow(epsilon0 * mc / 2, 0.5);
        kk *=
                (4 * PI() * pow(a, 4) * log(b)) /
                (pow(pow(nc / 2 * kb * c, 0.5) * (1 / e), 3) * pow(d, 7 / 2));
        return kk;
}
 
const long double
 
Ke = KoefK(z, De, Te, me),
Kc = KoefK(z, Di, Ti, mc),
Kh = KoefK(z, Dh, Th, mh), */
 
/* 23 */
//delta_ec = pow(Mve / Mvc, 2),
//delta_eh = pow(Mve / Mvh, 2),
//delta_ch = pow(Mvc / Mvh, 2);
 
/* 24 */
/*long double KoefH(
        const double & a,
        const double & b,
        const double & c)
{
        long double x;
        x = pow(a, 2) / (pow(b, 5) * pow(c, 3));
        return x;
}
 
const long double
 
He = KoefH(Mne, Mve, Mvc),
Hc = KoefH(Mnc, Mvc, Mve),
Hh = KoefH(Mnh, Mvh, Mve);*/
 
 
 
 
 
 
 
int main()
{
        setlocale(LC_ALL, "Russian");
        srand(time(0));
 
        /* 1 */
 
        long double
 
                ne = pow(10.0, 10),
                nc = pow(10.0, 16),
                nh = pow(10.0, 10),
 
                me = 9.11 * pow(10.0, -31),
                mc = 12.011 * 1.67 * pow(10.0, -27) - 9.11 * pow(10.0, -31),
                mh = 4.002 * 1.67 * pow(10.0, -27) - 9.11 * pow(10.0, -31),
 
                qe = -1.602176565 * pow(10.0, -19),
                qc = 1.602176565 * pow(10.0, -19),
                qh = 1.602176565 * pow(10.0, -19),
 
                //delta_t = pow(10, -12),
 
                rade = 2.81794 * pow(10.0, -15),
                radc = 31 * pow(10.0, -12),
                radh = 91 * pow(10.0, -12);
 
        /* 2 */
 
        long double
 
                Te = 4200,
                Ti = 4200,
                Th = 293,
 
                It = 150,
                //zc = 1,
 
                kb = 1.38 * pow(10.0, -23),
                shag_t = pow(10.0, -6),
                va = 0.44 * pow(10.0, -5),
                e = 1.602 * pow(10.0, -19),
 
                //z = 1, zi = 1, ze = 1, zh = 1,
 
                epsilon0 = 8.854187817 * pow(10.0, -12),
 
                gamma_e = (Te / Ti) / (me / mc),
                gamma_c = (Ti / Ti) / (mc / mc),
                gamma_h = (Th / Ti) / (mh / mc);
 
        /* 3 */
 
        long double
 
                xm = pow(10.0, -3),
                ym = 10 * pow(10.0, -3),
                zm = 10 * pow(10.0, -3),
 
                im = 5,
                km = 5,
                lm = 5,
 
                imm = 10,
                kmm = 10,
                lmm = 10,
 
                hx = xm / im,
                hy = ym / km,
                hz = zm / lm;
 
        /* 4 */
 
        long double
 
                hxx = xm / imm,
                hyy = ym / kmm,
                hzz = zm / lmm;
 
        /* 5 */
 
        long double
 
                Mc = nc * mc,
                Mh = nh * mh,
                Me = ne * me,
 
                rmax = xm,
                vmax = 10000,
                tmax = 1,
                NN = 300,
 
                Nr = NN,
                Nvr = NN,
                Nt = NN,
 
                hr = rmax / Nr;
 
        /* 6 */
 
        long double
 
                hvr = 2 * vmax / Nvr;
 
        /* 7 */
 
        long double
 
                ht = tmax / Nt,
                Tc = Ti,
                N = 1;
 
        /* 8 */
 
        long double
 
                Ml = pow((double)(kb * epsilon0 * Ti), 0.5) / pow((double)(e * nc), 0.5),
                Mnc = nc,
                Mne = ne,
                Mnh = nh,
 
                Mphi = kb * Ti / e;
 
        /* 9 */
 
        long double
 
                Mve = pow((double)(2 * kb * Te / me), 0.5),
                Mvc = pow((double)(2 * kb * Ti / mc), 0.5),
                Mvh = pow((double)(2 * kb * Th / mh), 0.5);
 
        /* 10 */
 
        long double
 
                Mee = Mphi / Ml;
 
        /* 11 */
 
        long double
 
                Mte = Ml / Mve,
                Mtc = Ml / Mvc,
                Mth = Ml / Mvh;
 
        /* 12 */
 
        long double
 
                Mfe = Mne / pow(Mve, 3),
                Mfc = Mnc / pow(Mvc, 3),
                Mfh = Mnh / pow(Mvh, 3);
 
        /* 13 */
 
        long double
 
                Mje = e * Mne * Mve;
 
        /* 14 */
 
        long double
 
                Mjc = e * Mnc * Mvc;
 
        /* 15 */
 
        long double
 
                Mjh = e * Mnh * Mvh;
 
        /* 16 */
 
        /* 17 */
 
        long double
 
                delta_t = 3.4 * 1 / (5.64 * pow(10.0, 4) * pow((double)(ne), 0.5));
 
        /* 18 */
 
        long double
 
                htn = ht / delta_t;
 
        /* 19 */
 
        long double
 
                hc = 1, // !!!!!
                hvrn = hvr / hc;
 
        /* 20 */
 
        long double
 
                hrn = hr / hc;
 
        /* 21 */
 
        //long double
 
                //zc = 1,
                //zh = 1,
                //ze = 1;
 
        /* 25 */
 
        //long double
 
                gamma_c = mc / me,
                gamma_e = me / me,
                gamma_h = mh / me;
 
        /* 26 */
        long double Vemax = pow((double)((2 * kb * Te) / me), 0.5);
 
        /* 27 */
        long double Vhmax = pow((double)((2 * kb * Th) / mh), 0.5);
 
        /* 28 */
        long double Vcmax = pow((double)((2 * kb * Tc) / mc), 0.5);
 
        /* 29 */
        long double Vhl = Vhmax - 0.1 * Vhmax;
 
        /* 30 */
        long double Vhpr = Vhmax + 0.1 * Vhmax;
 
        /* 31 */
        long double Vcl = Vcmax - 0.1 * Vcmax;
 
        /* 32 */
        long double Vcpr = Vcmax + 0.1 * Vcmax;
 
        /* 33 */
        long double Vel = Vemax - 0.1 * Vemax;
 
        /* 34 */
        long double Vepr = Vemax + 0.1 * Vemax;
 
        /* 35 */
        long double Kv2 = 0.005; // êîýôô. ñêîðîñòè
 
        /* 36 */
        long double Kv1 = pow((double)((3 - 2 * pow((double)(Kv2), 2))), 0.5);
 
        /* 37 */
 
        // Íà÷àëüíîå ïîëîæåíèå ÷àñòèö, øàã ÷àñòèöû
        
 
        long double        dd = rand() % 1;
		int
                je = 1,
                tt = 1;
 
        long double
                s1 = 0, s2 = 0, s3 = 0;
 //this 3000 stack overflow
        const int m = 1000, n =1000;
 
     long double
 
                **xe = new long double*[m],
                **ye = new long double*[m],
                **ze = new long double*[m],
 
                **xc = new long double*[m],
                **yc = new long double*[m],
                **zc = new long double*[m],
 
                **xh = new long double*[m],
                **yh = new long double*[m],
                **zh = new long double*[m],
                ////
                **rre = new long double*[m],
                **rrc = new long double*[m],
                **rrh = new long double*[m],
                ////
                **ve = new long double*[m],
                **veb = new long double*[m],
 
                **vxeb = new long double*[m],
                **vyeb = new long double*[m],
                **vzeb = new long double*[m],
                ////
                **xeb = new long double*[m],
                **yeb = new long double*[m],
                **zeb = new long double*[m],
 
                **xhb = new long double*[m],
                **yhb = new long double*[m],
                **zhb = new long double*[m],
 
                **xcb = new long double*[m],
                **ycb = new long double*[m],
                **zcb = new long double*[m],
                ////
                **rel = new long double*[m],
                **relb = new long double*[m],
 
                **rcl = new long double*[m],
                **rhl = new long double*[m],
                ///
                **rxel = new long double*[m],
                **ryel = new long double*[m],
                **rzel = new long double*[m],
 
                **rxhl = new long double*[m],
                **ryhl = new long double*[m],
                **rzhl = new long double*[m],
 
                ////
                **rxelb = new long double*[m],
                **ryelb = new long double*[m],
                **rzelb = new long double*[m],
 
                **rxclb = new long double*[m],
                **ryclb = new long double*[m],
                **rzclb = new long double*[m],
 
                **rxhlb = new long double*[m],
                **ryhlb = new long double*[m],
                **rzhlb = new long double*[m],
                ////
                **vxclb = new long double*[m],
                **vyclb = new long double*[m],
                **vzclb = new long double*[m],
 
                **vcl = new long double*[m],
                **vclb = new long double*[m];   
 
    for (int i = 0; i < n; i++)
    {
                xe[i] = new long double[n],
                ye[i] = new long double[n],
                ze[i] = new long double[n],
 
                xc[i] = new long double[n],
                yc[i] = new long double[n],
                zc[i] = new long double[n],
 
                xh[i] = new long double[n],
                yh[i] = new long double[n],
                zh[i] = new long double[n],
                ////
                rre[i] = new long double[n],
                rrc[i] = new long double[n],
                rrh[i] = new long double[n],
                ////
                ve[i] = new long double[n],
                veb[i] = new long double[n],
 
                vxeb[i] = new long double[n],
                vyeb[i] = new long double[n],
                vzeb[i] = new long double[n],
                ////
                xeb[i] = new long double[n],
                yeb[i] = new long double[n],
                zeb[i] = new long double[n],
 
                xhb[i] = new long double[n],
                yhb[i] = new long double[n],
                zhb[i] = new long double[n],
 
                xcb[i] = new long double[n],
                ycb[i] = new long double[n],
                zcb[i] = new long double[n],
                ////
                rel[i] = new long double[n],
                relb[i] = new long double[n],
 
                rcl[i] = new long double[n],
                rhl[i] = new long double[n],
                ///
                rxel[i] = new long double[n],
                ryel[i] = new long double[n],
                rzel[i] = new long double[n],
 
                rxhl[i] = new long double[n],
                ryhl[i] = new long double[n],
                rzhl[i] = new long double[n],
 
                ////
                rxelb[i] = new long double[n],
                ryelb[i] = new long double[n],
                rzelb[i] = new long double[n],
 
                rxclb[i] = new long double[n],
                ryclb[i] = new long double[n],
                rzclb[i] = new long double[n],
 
                rxhlb[i] = new long double[n],
                ryhlb[i] = new long double[n],
                rzhlb[i] = new long double[n],
                ////
                vxclb[i] = new long double[n],
                vyclb[i] = new long double[n],
                vzclb[i] = new long double[n],
 
                vcl[i] = new long double[n],
                vclb[i] = new long double[n];
    }
 
 
        HWND hwnd = GetConsoleWindow();
        HDC hdc = GetDC(hwnd); //ïîëó÷àåì DC(êîíòåêñò óñòðîéñòâà) äëÿ ðèñîâàíèÿ
        HPEN hpen1; //îáúÿâëÿåì îáúåêò ïåðî
        HGDIOBJ hpenOld, hbrushOld;
        HBRUSH hbrush; //îáúÿâëÿåì êèñòü
 
        hpen1 = CreatePen(PS_SOLID, 1, RGB(102, 0, 255)); //ëîãè÷åñêîå ïåðî ñ çàäàííûì ñòèëåì, øèðèíîé è öâåòîì
        hpenOld = (HPEN)SelectObject(hdc, hpen1);
 
        MoveToEx(hdc, 10, 10, NULL);//óñòàíàâëèâàåò òåêóùåé ïîçèöèåé óêàçàííóþ òî÷êó
        LineTo(hdc, 630, 10);
        LineTo(hdc, 630, 290);
        LineTo(hdc, 10, 290);
        LineTo(hdc, 10, 10);
 
        SetPixel(hdc, 43, 43, RGB(255, 0, 255));
        ReleaseDC(hwnd, hdc);
 
 
        for (int i = 1; i < imm; i++)
        {
                for (int k = 1; k < kmm; k++)
                {
                        for (int l = 1; l < lmm; l++)
                        {
                                s1 = hxx * (dd / pow(10.0, 12));
                                s2 = hyy * (dd / pow(10.0, 12));
                                s3 = hzz * (dd / pow(10.0, 12));
 
                                xe[je][tt] = (i - 1) * hxx + s1;
                                ye[je][tt] = (k - 1) * hyy + s2;
                                ze[je][tt] = (l - 1) * hzz + s3;
 
                                rre[je][tt] = pow(xe[je][tt], 2);
                                rre[je][tt] += pow(ye[je][tt], 2);
                                rre[je][tt] += pow(ze[je][tt], 2);
                                rre[je][tt] = pow((double)(rre[je][tt]), 0.5);
 
                                ve[je][tt] = pow((double)((3 / me) * kb * Te), 0.5);
 
                                veb[je][tt] = ve[je][tt] / Mve;
                                vxeb[je][tt] = (veb[je][tt] / 0.5) * Kv1;
                                vyeb[je][tt] = (veb[je][tt] / 0.5) * Kv2;
                                vzeb[je][tt] = (veb[je][tt] / 0.5) * Kv2;
 
                                xeb[je][tt] = xe[je][tt] / Ml;
                                yeb[je][tt] = ye[je][tt] / Ml;
                                zeb[je][tt] = ze[je][tt] / Ml;
 
                                rel[je][tt] = pow(xe[je][tt], 2);
                                rel[je][tt] += pow(ye[je][tt], 2);
                                rel[je][tt] += pow(ze[je][tt], 2);
                                rel[je][tt] = pow((double)(rel[je][tt]), 0.5);
 
                                relb[je][tt] = rel[je][tt] / Ml;
 
                                rxel[je][tt] = xe[je][tt];
                                ryel[je][tt] = ye[je][tt];
                                rzel[je][tt] = ze[je][tt];
 
 
                                rxelb[je][tt] = xeb[je][tt];
                                ryelb[je][tt] = yeb[je][tt];
                                rzelb[je][tt] = zeb[je][tt];
 
                                je++;
                        }
                }
        }
 
        int jme = je - 1;
 
        /* 38 */
        int jh = 1;
        tt = 1;
 
        for (int i = 1; i < imm; i++)
        {
                for (int k = 1; k < kmm; k++)
                {
                        for (int l = 1; l < lmm; l++)
                        {
                                s1 = hxx * (dd / pow(10.0, 12));
                                s2 = hyy * (dd / pow(10.0, 12));
                                s3 = hzz * (dd / pow(10.0, 12));
 
                                xh[jh][tt] = (i - 1) * hyy + s2;
                                yh[jh][tt] = (k - 1) * hyy + s2;
                                zh[jh][tt] = (l - 1) * hzz + s3;
 
                                rrh[jh][tt] = pow(xe[jh][tt], 2);
                                rrh[jh][tt] += pow(ye[jh][tt], 2);
                                rrh[jh][tt] += pow(ze[jh][tt], 2);
                                rrh[jh][tt] = pow((double)(rre[jh][tt]), 0.5);
 
                                rhl[jh][tt] = pow(xe[jh][tt], 2);
                                rhl[jh][tt] += pow(ye[jh][tt], 2);
                                rhl[jh][tt] += pow(ze[jh][tt], 2);
                                rhl[jh][tt] = pow((double)(rre[jh][tt]), 0.5);
 
                                vcl[jh][tt] = pow((double)((kb * Ti / mc)), 0.5);
                                vclb[jh][tt] = vcl[jh][tt] / Mvc;
 
                                xcb[jh][tt] = xh[jh][tt] / Ml;
                                ycb[jh][tt] = yh[jh][tt] / Ml;
                                zcb[jh][tt] = zh[jh][tt] / Ml;
 
                                rxhl[jh][tt] = xh[jh][tt];
                                ryhl[jh][tt] = yh[jh][tt];
                                rzhl[jh][tt] = zh[jh][tt];
 
                                rxhlb[jh][tt] = xhb[jh][tt];
                                ryhlb[jh][tt] = yhb[jh][tt];
                                rzhlb[jh][tt] = zhb[jh][tt];
 
                                jh++;
                        }
                }
        }
 
        /* 39 */
        int jc = 1;
        tt = 1;
 
        for (int i = 1; i < imm; i++)
        {
                for (int k = 1; k < kmm; k++)
                {
                        for (int l = 1; l < lmm; l++)
                        {
                                s1 = hxx * (dd / pow(10.0, 12));
                                s2 = hyy * (dd / pow(10.0, 12));
                                s3 = hzz * (dd / pow(10.0, 12));
    //this tt=1
                                xc[jc][tt] = 0;
                                yc[jc][tt] = (k - 1) * hyy + s2;
                                zc[jc][tt] = (l - 1) * hzz + s3;
 
                                rrc[jc][tt] = pow(xe[je][tt], 2);
                                rrc[jc][tt] += pow(ye[je][tt], 2);
                                rrc[jc][tt] += pow(ze[je][tt], 2);
                                rrc[jc][tt] = pow((double)(rre[je][tt]), 0.5);
 
                                vcl[jc][tt] = pow((double)((kb * Ti / mc)), 0.5);
                                vclb[jc][tt] = vcl[jc][tt] / Mvc;
 
                                xcb[jc][tt] = xc[jc][tt] / Ml;
                                ycb[jc][tt] = yc[jc][tt] / Ml;
                                zcb[jc][tt] = zc[jc][tt] / Ml;
 
                                rcl[jc][tt] = pow(xc[jc][tt], 2);
                                rcl[jc][tt] += pow(yc[jc][tt], 2);
                                rcl[jc][tt] += pow(zc[jc][tt], 2);
                                rcl[jc][tt] = pow((double)(rcl[jc][tt]), 0.5);
 
                                vxclb[jc][tt] = (vclb[jc][tt] / pow(3, 0.5) * Kv1);
                                vyclb[jc][tt] = (vclb[jc][tt] / pow(3, 0.5) * Kv2);
                                vzclb[jc][tt] = (vclb[jc][tt] / pow(3, 0.5) * Kv2);
 
                                rxclb[je][tt] = xcb[je][tt];
                                ryclb[je][tt] = ycb[je][tt];
                                rzclb[je][tt] = zcb[je][tt];
 
                                jc++;
                        }
                }
        }
 
        int jmc = jc - 1;
 
        _getch();
        return 0;
}
