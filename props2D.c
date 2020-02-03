/*
 *
 * Author: Mauricio Araya Polo
 * Date: 04/2015 - present
 *
 */

#include "props2D.h"

void bc2Dx( float *p, float *po, float *attr, float *left, float *right,
            int *prms, float *cfs, float *acfs, float *lcfs, float *rcfs )
{
    uint i, j, base, indxe, indx, nxe, jstencil;
    float tmp1, tmp2;

    TRANSLATE2D()
    
    nxe = nx + stencil + stencil;
    #pragma omp parallel for schedule(static) private(j, i, jstencil, base, tmp1, tmp2, indxe)
    #pragma ivdep
    for( j = stencil; j < nz + stencil; j++ )
    {
        base = j-stencil;

        /* left */
        tmp1 = attr[ base ] * dt;
        tmp1 *= tmp1;
        indxe = j + col;
        tmp2 = tmp1 * (acfs[0]* po[indxe      ] +
                       acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                       acfs[2]*(po[indxe - col] + po[indxe + col] ));
        p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                     - left[ base + nz ];
        indxe = j + col2;
        tmp2 = tmp1 * (acfs[0]* po[indxe      ] +
                       acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                       acfs[2]*(po[indxe - col] + po[indxe + col] ));
        p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                     - left[ base + 2*nz ];
        indxe = j + col3;
        tmp2 = tmp1 * (acfs[0]* po[indxe      ] +
                       acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                       acfs[2]*(po[indxe - col] + po[indxe + col] ));
        p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                     - left[ base + 3*nz ];
        p[ j ] = - ( lcfs[ base     ] * p[ j + col ] +
                     lcfs[ base + nz ] * p[ j + col2 ] +
                     lcfs[ base + 2*nz ] * po[ j ] +
                     lcfs[ base + 3*nz ] * po[ j + col ] +
                     lcfs[ base + 4*nz ] * po[ j + col2 ] +
                     lcfs[ base + 5*nz ] * left[ base     ] +
                     lcfs[ base + 6*nz ] * left[ base + nz ] +
                     lcfs[ base + 7*nz ] * left[ base + 2*nz ] ) ;

        /* right */
        tmp1 = attr[ base + (nx-1)*nz ] * dt;
        tmp1 *= tmp1;
	    i = nx + stencil;
            indxe = j + i*col;
            tmp2 = tmp1 * (acfs[0]* po[indxe      ] +
                           acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                           acfs[2]*(po[indxe - col] + po[indxe + col] ) );
            p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                         - right[ base + (nxe - i - 1)*nz ];
	    i++;
            indxe = j + i*col;
            tmp2 = tmp1 * (acfs[0]* po[indxe      ] +
                           acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                           acfs[2]*(po[indxe - col] + po[indxe + col] ) );
            p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                         - right[ base + (nxe - i - 1)*nz ];
	    i++;
            indxe = j + i*col;
            tmp2 = tmp1 * (acfs[0]* po[indxe      ] +
                           acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                           acfs[2]*(po[indxe - col] + po[indxe + col] ) );
            p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                         - right[ base + (nxe - i - 1)*nz ];
	    indxe = j + (nxe-1)*col;
        p[ indxe ] = - ( rcfs[ base     ] * p[ indxe - col ] +
                         rcfs[ base + nz ] * p[ indxe - col2 ] +
                         rcfs[ base + 2*nz ] * po[ indxe ] +
                         rcfs[ base + 3*nz ] * po[ indxe - col ] +
                         rcfs[ base + 4*nz ] * po[ indxe - col2 ] +
                         rcfs[ base + 5*nz ] * right[ base ] +
                         rcfs[ base + 6*nz ] * right[ base + nz ] +
                         rcfs[ base + 7*nz ] * right[ base + 2*nz] ) ;
    }
}

void fwd_step2D( float *p, float *po, float *attr, int *prms, float *cfs )
{
    float tmp2, tmp3[4000];
    uint indx, column, indxe, iz, ix;

    TRANSLATE2D()

    #pragma omp parallel for private(iz, ix, indx, indxe, tmp2, tmp3, column), schedule(runtime)
    for( ix = stencil; ix < nx + stencil; ix++ )
    {
        column = ix*col;
        // manual pre-fetching, non-contiguous memory access
        indxe = stencil + column;
        tmp3[0]  = c0 * po[ indxe ] +
                   cx1 * ( po[ indxe + col  ] + po[ indxe - col  ] ) +
                   cx2 * ( po[ indxe + col2 ] + po[ indxe - col2 ] ) +
                   cx3 * ( po[ indxe + col3 ] + po[ indxe - col3 ] ) +
                   cx4 * ( po[ indxe + col4 ] + po[ indxe - col4 ] ) ;
        #pragma ivdep
        for( iz = stencil + 1; iz < nz + stencil; iz++ )
        {
            indxe = iz + column;

            tmp3[iz-stencil] = c0 * po[ indxe ] +
                               cx1 * ( po[ indxe + col  ] + po[ indxe - col  ] ) +
                               cx2 * ( po[ indxe + col2 ] + po[ indxe - col2 ] ) +
                               cx3 * ( po[ indxe + col3 ] + po[ indxe - col3 ] ) +
                               cx4 * ( po[ indxe + col4 ] + po[ indxe - col4 ] ) ;
        }
        #pragma ivdep
        for( iz = stencil; iz < nz + stencil; iz++ )
        {
            indxe = iz + column;
            indx = (iz-stencil) + (ix-stencil) * nz;
            tmp2 = attr[ indx ] * dt;
            tmp2 *= tmp2;

            tmp3[iz-stencil] = tmp2*( tmp3[iz-stencil] +
                               cz1 * ( po[ indxe + 1 ] + po[ indxe - 1 ] ) +
                               cz2 * ( po[ indxe + 2 ] + po[ indxe - 2 ] ) +
                               cz3 * ( po[ indxe + 3 ] + po[ indxe - 3 ] ) +
                               cz4 * ( po[ indxe + 4 ] + po[ indxe - 4 ] )
                               );

            p[ indxe ] = po[ indxe ] + po[ indxe ] + tmp3[iz-stencil] - p[ indxe ];

        }
    }
}

void bc2Dz( float *p, float *po, float *attr, float *top, float *bottom,
            int *prms, float *cfs, float *acfs, float *tcfs, float *bcfs )
{
    uint k, j, column, indxe, indx, kstencil;
    float tmp, tmp2;

    TRANSLATE2D()

    #pragma omp parallel for private(k, j, kstencil, indx, indxe, tmp2, tmp, column), schedule(static)
    #pragma ivdep
    for( k = stencil; k < nx + stencil; k++ )
    {
        column = k*col;
        kstencil = k-stencil;
        /* top */
        indxe = 1 + column;
        indx = kstencil*nz;
        tmp = attr[ indx ] * dt;
        tmp *= tmp;
        tmp2 = tmp * ( acfs[0] * (po[ indxe ]) +
                       acfs[1] * (po[ indxe - 1   ] + po[ indxe + 1   ]) +
                       acfs[2] * (po[ indxe - col ] + po[ indxe + col ]) );
        p[ indxe ] = po[ indxe] + po[ indxe ] - top[ 1 + kstencil*stencil ] + tmp2;

        indxe = 2 + column;
        indx = 1 + kstencil*nz;
        tmp = attr[ indx ] * dt;
        tmp *= tmp;
        tmp2 = tmp * ( acfs[0] * (po[ indxe ]) +
                       acfs[1] * (po[ indxe - 1   ] + po[ indxe + 1   ]) +
                       acfs[2] * (po[ indxe - col ] + po[ indxe + col ]) );
        p[ indxe ] = po[ indxe] + po[ indxe ] - top[ 2 + kstencil*stencil ] + tmp2;

        indxe = 3 + column;
        indx = 2 + kstencil*nz;
        tmp = attr[ indx ] * dt;
        tmp *= tmp;
        tmp2 = tmp * ( acfs[0] * (po[ indxe ]) +
                       acfs[1] * (po[ indxe - 1   ] + po[ indxe + 1   ]) +
                       acfs[2] * (po[ indxe - col ] + po[ indxe + col ]) );
        p[ indxe ] = po[ indxe] + po[ indxe ] - top[ 3 + kstencil*stencil ] + tmp2;

        indxe = kstencil*8;
        indx = kstencil*stencil;
        p[ column ] = - ( tcfs[indxe    ] * p[ column + 1 ] +
                          tcfs[indxe + 1] * p[ column + 2 ] +
                          tcfs[indxe + 2] * po[ column ] +
                          tcfs[indxe + 3] * po[ column + 1] +
                          tcfs[indxe + 4] * po[ column + 2] +
                          tcfs[indxe + 5] * top[ indx ] +
                          tcfs[indxe + 6] * top[ indx + 1] +
                          tcfs[indxe + 7] * top[ indx + 2] );
        /* bottom */
        tmp = attr[ kstencil*nz + nz - 1 ] * dt;
        tmp *= tmp;
	    j = nz + stencil;
            indxe = column + j;
            tmp2 = tmp * ( acfs[0]* po[indxe      ] +
                           acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                           acfs[2]*(po[indxe - col] + po[indxe + col] ) );
            p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                         - bottom[ kstencil*stencil + col - j - 1 ];
	    j++;
            indxe = column + j;
            tmp2 = tmp * ( acfs[0]* po[indxe      ] +
                           acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                           acfs[2]*(po[indxe - col] + po[indxe + col] ) );
            p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                         - bottom[ kstencil*stencil + col - j - 1 ];
	    j++ ;
            indxe = column + j;
            tmp2 = tmp * ( acfs[0]* po[indxe      ] +
                           acfs[1]*(po[indxe - 1  ] + po[indxe + 1  ] ) +
                           acfs[2]*(po[indxe - col] + po[indxe + col] ) );
            p[ indxe ] = tmp2 + po[ indxe ] + po[ indxe ]
                         - bottom[ kstencil*stencil + col - j - 1 ];
        indxe = kstencil*8;
        indx = kstencil*stencil;
        p[ column + col - 1] = - ( bcfs[indxe    ] * p[ column + col - 2] +
                                   bcfs[indxe + 1] * p[ column + col - 3] +
                                   bcfs[indxe + 2] * po[ column + col - 1] +
                                   bcfs[indxe + 3] * po[ column + col - 2] +
                                   bcfs[indxe + 4] * po[ column + col - 3] +
                                   bcfs[indxe + 5] * bottom[ indx ] +
                                   bcfs[indxe + 6] * bottom[ indx + 1] +
                                   bcfs[indxe + 7] * bottom[ indx + 2] );
    }
}

