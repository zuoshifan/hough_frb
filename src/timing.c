#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define frb "FRB010125"
#define indir "../frb_data/"
#define file "5.4ms.cube"
#define Dl 600.0 /* lower bound of DM */
#define Dh 1000.0 /* uper bound of DM */
#define ND 2000 /* number of DM */
#define DM 4.5
#define dl Dl * DM
#define dh Dh * DM
#define dt 4.0
#define Nt 500
#define Nf 96


#define min(a,b)                                \
    ({ __typeof__ (a) _a = (a);                 \
        __typeof__ (b) _b = (b);                \
        _a > _b ? _b : _a; })

#define max(a,b)                                \
    ({ __typeof__ (a) _a = (a);                 \
        __typeof__ (b) _b = (b);                \
        _a > _b ? _a : _b; })

void insertion_sort(float *a, int n)
// Insertion sort array a of length n
{
    int i, j;
    float key;

    for(j = 1; j < n; j++)
    {
        key = a[j];
        // insert a[j] in the correct position a[0...(j-1)]
        i = j - 1;
        while ((i >= 0) && (a[i] > key))
        {
            a[i + 1] = a[i];
            i--;
        }
        a[i+1] = key;
    }
}


int partition(float *a, int n, float x)
// Partition array a of length n around x;
// Return the number of elements to the left of the pivot.
{
    int i, j;
    float tmp;

    // First find the pivot and place at the end
    for(i = 0; i < n; i++)
    {
        if(a[i] == x)
        /* if(fabs(a[i] - x) < 1.0e-8) */
        {
            a[i] = a[n-1];
            a[n-1] = x;
        }
    }

    i = 0;
    for(j = 0; j < (n-1); j++)
    {
        if(a[j] <= x)
        {
            tmp = a[j];
            a[j] = a[i];
            a[i] = tmp;
            i++;
        }
    }

    // Place the pivot in the correct position
    a[n-1] = a[i];
    a[i] = x;

    return i;
}

float ith_select(float *a, int i, int n)
// Select the ith element (indexed from 0) from the array of length n
// using the median of medians algorithm
// https://en.wikipedia.org/wiki/Median_of_medians
{
    int j, l;
    float tmp;

    if(n == 1)
    {
        return a[0];
    }

    int n_meds = 0;
    for(j = 0; j < n; j += 5)
    {
        l = min(5, n - j);
        insertion_sort(a + j, l);
        tmp = a[j/5];
        a[j/5] = a[j + l/2];
        a[j + l/2] = tmp;
        n_meds++;
    }

    float median_of_medians;
    if(n_meds > 1)
    {
        median_of_medians = ith_select(a, n_meds/2, n_meds);
    }
    else
    {
        median_of_medians = a[0];
    }

    int k = partition(a, n, median_of_medians);

    if(k == i)
    {
        return median_of_medians;
    }
    else if (i < k)
    {
        return ith_select(a, i, k);
    }
    else
    {
        return ith_select(a + k, i - k, n - k);
    }
}

/* float median(float *a, int n, int overwrite) */
/* // get the median of an array a of length n, */
/* // a will be changed if overwrite=1, else unchanged */
/* { */
/*     float *a1; */
/*     if(overwrite==0) /\* copy a to a1, so unchange a *\/ */
/*     { */
/*         a1 = (float *)malloc(n * sizeof(float)); */
/*         int i; */
/*         for(i=0; i<n; i++) */
/*             a1[i] = a[i]; */
/*     } */
/*     else /\* change a inplace *\/ */
/*     { */
/*         a1 = a; */
/*     } */

/*     return ith_select(a1, n/2, n); */
/* } */


int brute_dedisp(float **I_ft, float *time, float *freq, int nt, int nf, float dl1, float dh1, int nd, float **B)
// brute force dedispersion
{
    int i, di, fi, ti, nsh;
    float fl, fh, cf, tl, th, dt1, delta_t, dd, d;

    fl = freq[0];
    fh = freq[nf-1];
    cf = freq[nf/2]; /* central frequency */
    tl = time[0];
    th = time[nt-1];
    dt1 = (th - tl) / nt;

    float **I1 = (float **)malloc(nf * sizeof(float *));
    for (i=0; i<nf; i++)
         I1[i] = (float *)malloc(nt * sizeof(float));

    dd = (dh1 - dl1) / (nd - 1);
    for(di=0; di<nd; di++)
    {
        d = dl1 + dd * di;
        for(fi=0; fi<nf; fi++)
        {
            delta_t = d * (1.0/(freq[fi]*freq[fi]) - 1.0/(cf*cf));
            nsh = (int) delta_t/dt1;
            for(ti=0; ti<nt; ti++)
            {
                if(nsh>0) /* right shift */
                {
                    if(ti<nsh)
                    {
                        I1[fi][ti] = 0.0;
                    }
                    else
                    {
                        I1[fi][ti] = I_ft[fi][ti-nsh];
                    }
                }
                else if(nsh<0) /* left shitf */
                {
                    if(ti<nt+nsh)
                    {
                        I1[fi][ti] = I_ft[fi][ti-nsh];
                    }
                    else
                    {
                        I1[fi][ti] = 0.0;
                    }
                }
                else
                {
                    I1[fi][ti] = I_ft[fi][ti];
                }
            }
        }

        /* sum axis 0 of I1 to get B */
        for(ti=0; ti<nt; ti++)
            for(fi=0; fi<nf; fi++)
                B[di][ti] += I1[fi][ti];
    }

    /* free memory */
    for(i=0; i<nf; i++)
        free(I1[i]);
    free(I1);

    return 0;
}


int hough_transform(float **I_tf, float *time, float *freq, int nt, int nf, float dl1, float dh1, int nd, float t0l, float t0h, int nt0, float threshold, float **A)
// hough transform method
{
    int i, j, di, cnt=0, t0i;
    float dt0, d, dd;
    float *I, *Im, *t, *f;
    int len;
    float med, mad;

    dt0 = (t0h - t0l) / nt0;
    dd = (dh1 - dl1) / nd;

    len = nt * nf;
    /* compute median */
    I = (float *)malloc(len * sizeof(float));
    for(i=0; i<nt; i++)
        for(j=0; j<nf; j++)
            I[i*nf+j] = I_tf[i][j];
    med = ith_select(I, len/2, len);
    free(I);

    Im = (float *)malloc(len * sizeof(float));
    t = (float *)malloc(len * sizeof(float));
    f = (float *)malloc(len * sizeof(float));
    if(threshold > 0.0)
    {
        /* compute MAD */
        float *abs_diff = (float *)malloc(len * sizeof(float));
        float *abs_diff1 = (float *)malloc(len * sizeof(float));
        for(i=0; i<nt; i++)
            for(j=0; j<nf; j++)
            {
                abs_diff[i*nf+j] = fabs(I_tf[i][j] - med);
                abs_diff1[i*nf+j] = abs_diff[i*nf+j];
            }
        mad = ith_select(abs_diff, len/2, len); /* will change abs_diff */
        mad = mad / 0.6745;
        for(i=0; i<nt; i++)
            for(j=0; j<nf; j++)
            {
                if(abs_diff1[i] > threshold*mad)
                {
                    Im[cnt] = I_tf[i][j] - med; /* subtract mean */
                    t[cnt] = time[i];
                    f[cnt] = freq[j];
                    cnt++;
                }
            }
        free(abs_diff);
        free(abs_diff1);
    }
    else
    {
        for(i=0; i<nt; i++)
            for(j=0; j<nf; j++)
            {
                Im[cnt] = I_tf[i][j] - med; /* subtract mean */
                t[cnt] = time[i];
                f[cnt] = freq[j];
                cnt++;
            }
    }

    if (cnt == 0)
    {
        printf("No value that has threshold %f...\n", threshold);
    }
    else
    {
        for(i=0; i<cnt; i++)
        {
            for(di=0; di<nd; di++)
            {
                d = dl1 + dd*di;
                t0i = (int) ((-1.0/(f[i]*f[i]) + t[i]) - t0l) / dt0;
                t0i = max(0, t0i);
                t0i = min(nt0-1, t0i);
                A[t0i][di] += Im[i];
            }
        }
    }

    /* free memory */
    free(Im);
    free(t);
    free(f);

    return 0;
}


int main(int argc, char *argv[])
{
    int i, j;
    int fi, ti;
    float val;
    char fl_name[80];
    FILE *fp;

    float **I_ft = (float **)malloc(Nf * sizeof(float *));
    for (i=0; i<Nf; i++)
         I_ft[i] = (float *)malloc(Nt * sizeof(float));

    float **I_tf = (float **)malloc(Nt * sizeof(float *));
    for (i=0; i<Nt; i++)
         I_tf[i] = (float *)malloc(Nf * sizeof(float));
    float *t = (float *)malloc(Nt * sizeof(float));
    float *f = (float *)malloc(Nf * sizeof(float));

    sprintf(fl_name, "%s%s/%s", indir, frb, file);
    printf("%s\n", fl_name);
    fp = fopen(fl_name, "r");
    if(fp == NULL)
    {
        printf("File %s not found!", fl_name);
    }
    /* read in data */
    for (i=0; i<Nf; i++)
        for (j=0; j<Nt; j++)
        {
            fscanf(fp, "%d %d %f", &fi, &ti, &val);
            /* printf("%d, %d, %f\n", fi, ti, val); */
            I_ft[fi][ti] = val;
            I_tf[ti][fi] = val;
        }
    fclose(fp);

    sprintf(fl_name, "%s%s/%s", indir, frb, "freq.txt");
    printf("%s\n", fl_name);
    fp = fopen(fl_name, "r");
    if(fp == NULL)
    {
        printf("File %s not found!", fl_name);
    }
    /* read in freq */
    i = 0;
    while(i<Nf && fscanf(fp, "%f,", &val)!=EOF)
        f[Nf - (++i)] = 0.001*val; /* GHz, ascending */
    fclose(fp);

    /* generate time array */
    for (i=0; i<Nt; i++)
        t[i] = dt * i; /* ms */

    /* for timing */
    int tm_use, num=10;
    float tms[7]; /* to record all timings, in us */
    struct timeval start, end;

    /* brute force dedispersion */
    float **B = (float **)malloc(ND * sizeof(float *));
    for (i=0; i<ND; i++)
         B[i] = (float *)malloc(Nt * sizeof(float));

    tm_use = 0;
    for(i=0; i<num; i++)
    {
        gettimeofday( &start, NULL );
        brute_dedisp(I_ft, t, f, Nt, Nf, dl, dh, ND, B);
        gettimeofday( &end, NULL );
        tm_use += 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    }
    printf("time: %f us\n", ((float) tm_use)/num);
    tms[0] = ((float) tm_use)/num;
    /* free B */
    for(i=0; i<ND; i++)
        free(B[i]);
    free(B);

    /* hough transform dedispersion */
    float fl = f[0], fh = f[Nf-1];
    float tl = f[0], th = t[Nt-1];
    float t0l = -1.0/(fl*fl) * dh + tl;
    float t0h = -1.0/(fh*fh) * dl + th;
    int nt0 = Nt;
    float **A = (float **)malloc(nt0 * sizeof(float *));
    for (i=0; i<nt0; i++)
         A[i] = (float *)malloc(ND * sizeof(float));

    float threshold;
    float thresholds[11] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
    int tdi = 0;
    for (tdi=0; tdi<11; tdi++)
    {
        /* initialize A to zero */
        for(i=0; i<nt0; i++)
            for(j=0; j<ND; j++)
                A[i][j] = 0.0;

        threshold = thresholds[tdi];
        tm_use = 0;
        for(i=0; i<num; i++)
        {
            gettimeofday( &start, NULL );
            hough_transform(I_tf, t, f, Nt, Nf, dl, dh, ND, t0l, t0h, nt0, threshold, A);
            gettimeofday( &end, NULL );
            tm_use += 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        }
        printf("threshold %f, time: %f us\n", threshold, ((float) tm_use)/num);
        tms[tdi+1] = ((float) tm_use)/num;
    }
    /* free A */
    for(i=0; i<nt0; i++)
        free(A[i]);
    free(A);

    /* print all timings */
    printf("All timings: ");
    for(i=0; i<12; i++)
        printf("%f ", tms[i]);
    printf("\n");


    /* free memory */
    for(i=0; i<Nf; i++)
        free(I_ft[i]);
    free(I_ft);
    for(i=0; i<Nt; i++)
        free(I_tf[i]);
    free(I_tf);
    free(t);
    free(f);

    return 0;
}
