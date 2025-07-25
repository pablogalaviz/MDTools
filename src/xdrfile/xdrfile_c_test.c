/* Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* Get HAVE_RPC_XDR_H, F77_FUNC from config.h if available */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD
#include "unistd.h"
#endif

/* define M_PI for C89 */
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/* get fixed-width types if we are using ANSI C99 */
#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif (defined HAVE_INTTYPES_H)
#include <inttypes.h>
#endif

#ifdef HAVE_RPC_XDR_H
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#endif

#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"
#include <float.h>
#include <time.h>

/* This program tests reading and writing to XDR files */

static void _die(char* msg, int line, char* file) {
    fprintf(stderr, "Fatal error: %s\n", msg);
    fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
    exit(1);
}
#define die(msg) _die(msg, __LINE__, __FILE__)

static void _die_r(char* msg, int result, int line, char* file) {
    fprintf(stderr, "Fatal error: %s\n", msg);
    fprintf(stderr, "result = %d\n", result);
    fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
    exit(1);
}
#define die_r(msg, res) _die_r(msg, res, __LINE__, __FILE__)

static void test_xtc() {
    char* testfn = "test.xtc";
    XDRFILE* xd;
    int result, i, j, k, nframes = 13;
    int natoms2, natoms1 = 173;
    int step2, step1 = 1993;
    float time2, time1 = 1097.23;
    matrix box2, box1;
    rvec *x2, *x1;
    float prec2, prec1 = 1000;
    float toler = 1e-3;

    printf("Testing xtc functionality:");
    for (i = 0; (i < DIM); i++) {
        for (j = 0; (j < DIM); j++) {
            box1[i][j] = (i + 1) * 3.7 + (j + 1);
        }
    }
    x1 = calloc(natoms1, sizeof(*x1));
    if (NULL == x1) {
        die("Allocating memory for x1 in test_xtc");
    }

    for (i = 0; (i < natoms1); i++) {
        for (j = 0; (j < DIM); j++) {
            x1[i][j] = (i + 1) * 3.7 + (j + 1);
        }
    }
    xd = xdrfile_open(testfn, "w");
    if (NULL == xd) {
        die("Opening xdrfile for writing");
    }
    for (k = 0; (k < nframes); k++) {
        result = write_xtc(xd, natoms1, step1 + k, time1 + k, box1, x1, prec1);
        if (0 != result) {
            die_r("Writing xtc file", result);
        }
    }
    xdrfile_close(xd);

    result = read_xtc_natoms(testfn, &natoms2);
    if (exdrOK != result) {
        die_r("read_xtc_natoms", result);
    }
    if (natoms2 != natoms1) {
        die("Number of atoms incorrect when reading trr");
    }
    x2 = calloc(natoms2, sizeof(x2[0]));
    if (NULL == x2) {
        die("Allocating memory for x2");
    }

    xd = xdrfile_open(testfn, "r");
    if (NULL == xd) {
        die("Opening xdrfile for reading");
    }

    k = 0;
    do {
        result = read_xtc(xd, natoms2, &step2, &time2, box2, x2, &prec2);
        if (exdrENDOFFILE != result) {
            if (exdrOK != result) {
                die_r("read_xtc", result);
            }
            if (natoms2 != natoms1) {
                die("natoms2 != natoms1");
            }
            if (step2 - step1 != k) {
                die("incorrect step");
            }
            if (fabs(time2 - time1 - k) > toler) {
                die("incorrect time");
            }
            if (fabs(prec2 - prec1) > toler) {
                die("incorrect precision");
            }
            for (i = 0; (i < DIM); i++) {
                for (j = 0; (j < DIM); j++) {
                    if (fabs(box2[i][j] - box1[i][j]) > toler) {
                        die("box incorrect");
                    }
                }
            }
            for (i = 0; (i < natoms1); i++) {
                for (j = 0; (j < DIM); j++) {
                    if (fabs(x2[i][j] - x1[i][j]) > toler) {
                        die("x incorrect");
                    }
                }
            }
        }
        k++;
    } while (result == exdrOK);

    free(x1);
    free(x2);

    xdrfile_close(xd);
#ifdef HAVE_UNISTD
    unlink(testfn);
#endif
    printf(" PASSED\n");
}

static void test_trr() {
    char* testfn = "test.trr";
    XDRFILE* xd;
    int result, i, j, k, nframes = 13;
    int natoms2, natoms1 = 173;
    int step2, step1 = 1993;
    float time2, time1 = 1097.23;
    matrix box2, box1;
    rvec *x2, *x1;
    float lambda2, lambda1 = 0.4;
    float toler = 1e-3;
    uint8_t flag_buf = 0;

    printf("Testing trr functionality:");
    for (i = 0; (i < DIM); i++) {
        for (j = 0; (j < DIM); j++) {
            box1[i][j] = (i + 1) * 3.7 + (j + 1);
        }
    }
    x1 = calloc(natoms1, sizeof(*x1));
    if (NULL == x1) {
        die("Allocating memory for x1 in test_xtc");
    }

    for (i = 0; (i < natoms1); i++) {
        for (j = 0; (j < DIM); j++) {
            x1[i][j] = (i + 1) * 3.7 + (j + 1);
        }
    }
    xd = xdrfile_open(testfn, "w");
    if (NULL == xd) {
        die("Opening trr file for writing");
    }
    for (k = 0; (k < nframes); k++) {
        result = write_trr(xd, natoms1, step1 + k, time1 + k, lambda1, box1, x1, NULL, NULL);
        if (0 != result) {
            die_r("Writing trr file", result);
        }
    }
    xdrfile_close(xd);

    result = read_trr_natoms(testfn, &natoms2);
    if (exdrOK != result) {
        die_r("read_trr_natoms", result);
    }
    if (natoms2 != natoms1) {
        die("Number of atoms incorrect when reading trr");
    }
    x2 = calloc(natoms2, sizeof(x2[0]));
    if (NULL == x2) {
        die("Allocating memory for x2");
    }

    xd = xdrfile_open(testfn, "r");
    if (NULL == xd) {
        die("Opening trr file for reading");
    }

    for (k = 0; (k < nframes); k++) {
        result = read_trr(xd, natoms2, &step2, &time2, &lambda2, box2, x2, NULL, NULL, &flag_buf);
        if (exdrOK != result) {
            die_r("read_xtc", result);
        }
        if (natoms2 != natoms1) {
            die("natoms2 != natoms1");
        }
        if (step2 - step1 != k) {
            die("incorrect step");
        }
        if (fabs(time2 - time1 - k) > toler) {
            die("incorrect time");
        }
        if (fabs(lambda2 - lambda2) > toler) {
            die("incorrect lambda");
        }
        for (i = 0; (i < DIM); i++) {
            for (j = 0; (j < DIM); j++) {
                if (fabs(box2[i][j] - box1[i][j]) > toler) {
                    die("box incorrect");
                }
            }
        }
        for (i = 0; (i < natoms1); i++) {
            for (j = 0; (j < DIM); j++) {
                if (fabs(x2[i][j] - x1[i][j]) > toler) {
                    die("x incorrect");
                }
            }
        }
    }

    free(x1);
    free(x2);

    xdrfile_close(xd);
#ifdef HAVE_UNISTD
    unlink(testfn);
#endif
    printf(" PASSED\n");
}

static void test_basic() {
    float test_ii; /* 7 significant digits */
#define EPSILON_1 1e-7
#define EPSILON_2 1e-4

    printf("Testing basic xdrfile library:");
    for (test_ii = 1.0e1; test_ii < 1.0e2; (test_ii = test_ii + pow(M_PI, 0.00011))) {

#define BUFLEN 37
        XDRFILE* xfp;
        int i, j, k, len, ncoord = BUFLEN / 3;
        char ptr[BUFLEN], *buf = "abcdefghijklmnopqrstuvwxyz";
        char* testfn = "test.xdr";
        unsigned char uptr[BUFLEN];
        short sptr[BUFLEN], sptr2[BUFLEN];
        unsigned short usptr[BUFLEN], usptr2[BUFLEN];
        int iptr[BUFLEN], iptr2[BUFLEN];
        unsigned int uiptr[BUFLEN], uiptr2[BUFLEN];
        float fptr[BUFLEN], fptr2[BUFLEN];
        double dptr[BUFLEN], dptr2[BUFLEN];
        char optr[BUFLEN], optr2[BUFLEN];
#define NPREC 1

        float fprec[] = {234.45};
        double dprec[] = {234.45};

        /* Can not write a string that's on the stack since all data is
           treated as variables.
         */
        len = strlen(buf) + 1;
        if (len >= BUFLEN) {
            die("Increase BUFLEN");
        }
        strcpy(ptr, buf);
        strcpy((char*)uptr, buf);
        /* Initiate arrays */
        for (i = 0; (i < BUFLEN); i++) {
            fptr[i] = cos(i * 13.0 / M_PI);
            dptr[i] = sin(i * 13.0 / M_PI);
            iptr[i] = (int)(floor(dptr[i] * 1000));
            uiptr[i] = (unsigned int)(floor(dptr[i] * 1000) + 1001);
        }
        /* Initiate opaque array */
        memcpy(optr, dptr, BUFLEN);

        /*************************************/
        /*           WRITING BIT             */
        /*************************************/

        if ((xfp = xdrfile_open("test.xdr", "w")) == NULL) {
            die("Can not open file for writing");
        }

        if (xdrfile_write_char(ptr, len, xfp) != len) {
            die("Writing char string");
        }
        if (xdrfile_write_uchar(uptr, len, xfp) != len) {
            die("Writing uchar string");
        }
        if (xdrfile_write_short(sptr, BUFLEN, xfp) != BUFLEN) {
            die("Writing short array");
        }
        if (xdrfile_write_ushort(usptr, BUFLEN, xfp) != BUFLEN) {
            die("Writing ushort array");
        }
        if (xdrfile_write_int(iptr, BUFLEN, xfp) != BUFLEN) {
            die("Writing int array");
        }
        if (xdrfile_write_uint(uiptr, BUFLEN, xfp) != BUFLEN) {
            die("Writing uint array");
        }
        if (xdrfile_write_float(fptr, BUFLEN, xfp) != BUFLEN) {
            die("Writing float array");
        }
        if (xdrfile_write_double(dptr, BUFLEN, xfp) != BUFLEN) {
            die("Writing double array");
        }
        if (xdrfile_write_string(buf, xfp) != len) {
            die("Writing string");
        }
        if (xdrfile_write_opaque(optr, BUFLEN, xfp) != BUFLEN) {
            die("Writing opaque");
        }
        for (k = 0; (k < NPREC); k++) {
            if (xdrfile_compress_coord_float(fptr, ncoord, fprec[k], xfp) != ncoord) {
                die("Writing compress_coord_float");
            }
            if (xdrfile_compress_coord_double(dptr, ncoord, dprec[k], xfp) != ncoord) {
                die("Writing compress_coord_double");
            }
        }
        if (xdrfile_close(xfp) != 0) {
            die("Can not close xdr file");
        }

        /*************************************/
        /*          READING BIT              */
        /*************************************/
        if ((xfp = xdrfile_open(testfn, "r")) == NULL) {
            die("Can not open file for reading");
        }

        if ((xdrfile_read_char(ptr, len, xfp)) != len) {
            die("Not the right number of chars read from string");
        }
        if (strcmp(ptr, buf) != 0) {
            printf("did not read the expected chars");
        }
        if (xdrfile_read_uchar(uptr, len, xfp) != len) {
            die("Not the right number of uchars read from string");
        }
        if (strcmp((char*)uptr, buf) != 0) {
            printf("did not read the expected uchars");
        }
        if (xdrfile_read_short(sptr2, BUFLEN, xfp) != BUFLEN) {
            die("Reading short array");
        }
        for (i = 0; (i < BUFLEN); i++) {
            if (sptr2[i] != sptr[i]) {
                fprintf(stderr, "i: %5d, wrote: %10d, read: %10d\n", i, sptr[i], sptr2[i]);
                die("Comparing short array");
            }
        }
        if (xdrfile_read_ushort(usptr2, BUFLEN, xfp) != BUFLEN) {
            die("Reading ushort array");
        }
        for (i = 0; (i < BUFLEN); i++) {
            if (usptr2[i] != usptr[i]) {
                fprintf(stderr, "i: %5d, wrote: %10d, read: %10d\n", i, usptr[i], usptr2[i]);
                die("Comparing ushort array");
            }
        }
        if (xdrfile_read_int(iptr2, BUFLEN, xfp) != BUFLEN) {
            die("Reading int array");
        }
        for (i = 0; (i < BUFLEN); i++) {
            if (iptr2[i] != iptr[i]) {
                fprintf(stderr, "i: %5d, wrote: %10d, read: %10d\n", i, iptr[i], iptr2[i]);
                die("Comparing int array");
            }
        }
        if (xdrfile_read_uint(uiptr2, BUFLEN, xfp) != BUFLEN) {
            die("Reading uint array");
        }
        for (i = 0; (i < BUFLEN); i++) {
            if (uiptr2[i] != uiptr[i]) {
                fprintf(stderr, "i: %5d, wrote: %10d, read: %10d\n", i, uiptr[i], uiptr2[i]);
                die("Comparing uint array");
            }
        }
        if (xdrfile_read_float(fptr2, BUFLEN, xfp) != BUFLEN) {
            die("Reading float array");
        }
        for (i = 0; (i < BUFLEN); i++) {
            if (fptr2[i] != fptr[i]) {
                fprintf(stderr, "i: %5d, wrote: %12g, read: %12g\n", i, fptr[i], fptr2[i]);
                die("Comparing float array");
            }
        }
        if (xdrfile_read_double(dptr2, BUFLEN, xfp) != BUFLEN) {
            die("Reading double array");
        }
        for (i = 0; (i < BUFLEN); i++) {
            if (dptr2[i] != dptr[i]) {
                fprintf(stderr, "i: %5d, wrote: %12g, read: %12g\n", i, dptr[i], dptr2[i]);
                die("Comparing double array");
            }
        }
        if (xdrfile_read_string(ptr, BUFLEN, xfp) != len) {
            die("Reading string");
        }
        if (strcmp(ptr, buf) != 0) {
            die("Comparing strings");
        }
        if (xdrfile_read_opaque(optr2, BUFLEN, xfp) != BUFLEN) {
            die("Reading opaque array");
        }
        for (i = 0; (i < BUFLEN); i++) {
            if (optr2[i] != optr[i]) {
                fprintf(stderr, "i: %5d, wrote: %2d, read: %2d\n", i, optr[i], optr2[i]);
                die("Comparing opaque array");
            }
        }
        for (k = 0; (k < NPREC); k++) {
            float ff, fx;
            double dd, dx;
            int nc = ncoord;
            if (xdrfile_decompress_coord_float(fptr2, &nc, &ff, xfp) != ncoord) {
                die("Reading compress_coord_float");
            }
            if (fabs(ff - fprec[k]) > EPSILON_1) {
                printf("Found precision %f, expected %f\n", ff, fprec[k]);
                die("Float precision");
            }
            if (ff <= 0) {
                ff = 1000;
            }

            for (i = 0; (i < ncoord); i++) {
                for (j = 0; (j < 3); j++) {
                    fx = rint(fptr[3 * i + j] * ff) / ff;
                    if (fabs(fx - fptr2[3 * i + j]) > EPSILON_1) {
                        printf("prec: %10g, i: %3d, j: %d, fx: %10g, fptr2: %12g, fptr: %12g\n", ff,
                               i, j, fx, fptr2[3 * i + j], fptr[3 * i + j]);
                        die("Reading decompressed float coordinates");
                    }
                }
            }
            if (xdrfile_decompress_coord_double(dptr2, &nc, &dd, xfp) != ncoord) {
                die("Reading compress_coord_double");
            }

            if (fabs(dd - dprec[k]) > EPSILON_2) {
                die("Double precision");
            }

            for (i = 0; (i < ncoord); i++) {
                for (j = 0; (j < 3); j++) {
                    dx = rint(dptr[3 * i + j] * dd) / dd;
                    if (fabs(dx - dptr2[3 * i + j]) > EPSILON_2) {
                        printf("prec: %10g, i: %3d, j: %d, dx: %10g, dptr2: %12g, dptr: %12g\n", dd,
                               i, j, dx, dptr2[3 * i + j], dptr[3 * i + j]);
                        die("Reading decompressed double coordinates");
                    }
                }
            }
        }

        if (xdrfile_close(xfp) != 0) {
            die("Can not close xdr file");
        }
    }
#ifdef HAVE_UNISTD
    unlink(testfn);
#endif
    printf(" PASSED\n");
}

static void test_trr_offsets() {
    const char* testfn = "../test_data/traj_n3.trr";
    int result = 0;
    int natoms = 0;
    unsigned long i, nframes = 0;
    int64_t* offsets = NULL;
    const int natoms_ref = 3;
    const unsigned long nframes_ref = 101;
    const int64_t framesize_ref = 156;

    printf("Testing trr offset functionality: ");

    result = read_trr_header(testfn, &natoms, &nframes, &offsets);
    if (exdrOK != result) {
        die_r("read_trr_header", result);
    }

    if (natoms != natoms_ref) {
        printf("Found number of atoms %d, expected %d\n", natoms, natoms_ref);
        die("Number of atoms incorrect when reading trr");
    }

    if (nframes != nframes_ref) {
        printf("Found number of frames %lu, expected %lu\n", nframes, nframes_ref);
        die("Number of frames incorrect when reading trr");
    }

    for (i = 0; i < nframes; i++) {
        if (offsets[i] != i * framesize_ref) {
            printf("Found offsets[%lu] %ld, expected %lu\n", i, offsets[i], i * framesize_ref);
            die("Offset incorrect when reading trr");
        }
    }

    free(offsets);

#ifdef HAVE_UNISTD
    unlink(testfn);
#endif
    printf(" PASSED\n");
}

static void test_xtc_uncompressed_offsets() {
    const char* testfn = "../test_data/traj_comp_n3.xtc";
    int result = 0;
    int natoms = 0;
    unsigned long i, nframes = 0;
    int64_t* offsets = NULL;
    const int natoms_ref = 3;
    const unsigned long nframes_ref = 101;
    const int64_t framesize_ref = 92;

    printf("Testing uncompressed xtc offset functionality: ");

    result = read_xtc_header(testfn, &natoms, &nframes, &offsets);
    if (exdrOK != result) {
        die_r("read_xtc_header", result);
    }

    if (natoms != natoms_ref) {
        printf("Found number of atoms %d, expected %d\n", natoms, natoms_ref);
        die("Number of atoms incorrect when reading xtc");
    }

    if (nframes != nframes_ref) {
        printf("Found number of frames %lu, expected %lu\n", nframes, nframes_ref);
        die("Number of frames incorrect when reading xtc");
    }

    for (i = 0; i < nframes; i++) {
        if (offsets[i] != i * framesize_ref) {
            printf("Found offsets[%lu] %ld, expected %lu\n", i, offsets[i], i * framesize_ref);
            die("Offset incorrect when reading xtc");
        }
    }

    free(offsets);

#ifdef HAVE_UNISTD
    unlink(testfn);
#endif
    printf(" PASSED\n");
}

static void test_xtc_compressed_offsets() {
    const char* testfn = "../test_data/traj_comp_n11.xtc";
    XDRFILE* xd;
    int magic = 0;
    int result = 0;
    int natoms = 0;
    unsigned long i, nframes = 0;
    int64_t* offsets = NULL;
    const int natoms_ref = 11;
    const unsigned long nframes_ref = 101;
    const int magic_ref = 1995;

    printf("Testing compressed xtc offset functionality: ");

    result = read_xtc_header(testfn, &natoms, &nframes, &offsets);
    if (exdrOK != result) {
        die_r("read_xtc_header", result);
    }

    if (natoms != natoms_ref) {
        printf("Found number of atoms %d, expected %d\n", natoms, natoms_ref);
        die("Number of atoms incorrect when reading xtc");
    }

    if (nframes != nframes_ref) {
        printf("Found number of frames %lu, expected %lu\n", nframes, nframes_ref);
        die("Number of frames incorrect when reading xtc");
    }

    xd = xdrfile_open(testfn, "r");
    if (NULL == xd) {
        die("Opening xdrfile for reading");
    }

    for (i = 0; i < nframes; i++) {
        xdr_seek(xd, offsets[i], SEEK_SET);
        if (xdrfile_read_int(&magic, 1, xd) != 1) {
            die("xdrfile_read_int");
        }
        if (magic != magic_ref) {
            printf("Found magic %d, expected %d\n", magic, magic_ref);
            die("Incorrect magic when reading xtc");
        }
    }

    xdrfile_close(xd);

    free(offsets);

#ifdef HAVE_UNISTD
    unlink(testfn);
#endif
    printf(" PASSED\n");
}

static void test_trr_flag() {
    const char* testfn = "../test_data/traj_n3_xvf.trr";
    XDRFILE* xd;
    int result = 0;
    int natoms = 0;
    unsigned long i, nframes = 0;
    int64_t* offsets = NULL;
    int step = 0;
    float time = 0;
    float lambda = 0;
    matrix box;
    const int natoms_ref = 3;
    const unsigned long nframes_ref = 75;
    uint8_t flag = 0;

    printf("Testing trr flag functionality: ");

    result = read_trr_header(testfn, &natoms, &nframes, &offsets);
    if (exdrOK != result) {
        die_r("read_trr_header", result);
    }

    if (natoms != natoms_ref) {
        printf("Found number of atoms %d, expected %d\n", natoms, natoms_ref);
        die("Number of atoms incorrect when reading trr");
    }

    if (nframes != nframes_ref) {
        printf("Found number of frames %lu, expected %lu\n", nframes, nframes_ref);
        die("Number of frames incorrect when reading trr");
    }

    xd = xdrfile_open(testfn, "r");
    if (NULL == xd) {
        die("Opening trr file for reading");
    }

    for (i = 0; i < nframes; i++) {
        result = read_trr(xd, natoms, &step, &time, &lambda, box, NULL, NULL, NULL, &flag);
        if (exdrOK != result) {
            die_r("read_trr", result);
        }

        if (!(flag & TRR_HAS_BOX)) {
            printf("Expected box in frame %d\n", i);
            die("Box flag incorrect");
        }
        if (step % 2000 == 0 && !(flag & TRR_HAS_POSITIONS)) {
            printf("Expected positions in frame %d\n", i);
            die("Position flag incorrect");
        }
        if (step % 3000 == 0 && !(flag & TRR_HAS_VELOCITIES)) {
            printf("Expected velocities in frame %d\n", i);
            die("Velocity flag incorrect");
        }
        if (step % 5000 == 0 && !(flag & TRR_HAS_FORCES)) {
            printf("Expected forces in frame %d\n", i);
            die("Force flag incorrect");
        }
    }

    free(offsets);

#ifdef HAVE_UNISTD
    unlink(testfn);
#endif
    printf(" PASSED\n");
}

int main(int argc, char* argv[]) {
    /* Test basic stuff */
    test_basic();
    /* Now test writing a complete xtc file */
    test_xtc();

    test_trr();

    /* Test offsets of trr file by comparing to hard coded framesize */
    test_trr_offsets();

    /* Test offsets of uncompressed xtc file by comparing to hard coded framesize */
    test_xtc_uncompressed_offsets();

    /* Test offsets of compressed xtc file by using seek*/
    test_xtc_compressed_offsets();

    /* Test flag for data fields (box, x, v, f) in trr file */
    test_trr_flag();

    return 0;
}
