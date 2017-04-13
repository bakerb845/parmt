#include <stdio.h>
#include <stdlib.h>
#include "parmt_utils.h"
#ifdef PARMT_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#include "iscl/array/array.h"

/*!
 * @brief The iloc'th location and iobs'th observation this sets the forward
 *        modeling matrix G s.t. Gm = u
 *
 * @param[in] iobs    observation number
 * @param[in] iloc    location number
 * @param[in] npmax   leading dimension of G
 * @param[in] data    data data structure with Green's functions which will be
 *                    set on G
 *
 * @param[out] G      column major modeling matrix s.t. Gm=u where
 *                    G is packed
 *                     \f$ \{ G_{xx}, G_{yy}, G_{zz},
 *                            G_{xy}, G_{xz}, G_{yz} \} \f$
 *                    and m is in NED ordering [npmax x 6]
 *
 * @result 0 indicates success
 *
 */
int parmt_utils_setDataOnG(const int iobs, const int iloc, const int npmax,
                           const struct parmtData_struct data,
                           double *__restrict__ G)
{
    const char *fcnm = "parmt_utils_setDataOnG\0";
    int k;
    if (iobs < 0 || iobs >= data.nobs)
    {   
        printf("%s: Error invalid observation\n", fcnm);
        return -1; 
    }   
    if (iloc < 0 || iloc >= data.nlocs)
    {   
        printf("%s: Error invalid location\n", fcnm);
        return -1; 
    }   
    k = iobs*data.nlocs + iloc;
    if (npmax > data.sacGxx[k].npts || npmax > data.sacGyy[k].npts ||
        npmax > data.sacGzz[k].npts || npmax > data.sacGxy[k].npts ||
        npmax > data.sacGxz[k].npts || npmax > data.sacGyz[k].npts)
    {   
        printf("%s: npmax is too small\n", fcnm);
        return -1; 
    }   
    if (data.sacGxx[k].data == NULL || data.sacGyy[k].data == NULL ||
        data.sacGzz[k].data == NULL || data.sacGxy[k].data == NULL ||
        data.sacGxz[k].data == NULL || data.sacGyz[k].data == NULL)
    {   
        printf("%s: Greens function is null\n", fcnm);
        return -1; 
    }   
    // Get the Green's functions onto the matrix G
    array_zeros64f_work(6*npmax, G); 
    cblas_dcopy(data.sacGxx[k].header.npts,
                data.sacGxx[k].data, 1, &G[0*npmax], 1); 
    cblas_dcopy(data.sacGyy[k].header.npts,
                data.sacGyy[k].data, 1, &G[1*npmax], 1); 
    cblas_dcopy(data.sacGzz[k].header.npts,
                data.sacGzz[k].data, 1, &G[2*npmax], 1); 
    cblas_dcopy(data.sacGxy[k].header.npts,
                data.sacGxy[k].data, 1, &G[3*npmax], 1); 
    cblas_dcopy(data.sacGxz[k].header.npts,
                data.sacGxz[k].data, 1, &G[4*npmax], 1); 
    cblas_dcopy(data.sacGyz[k].header.npts,
                data.sacGyz[k].data, 1, &G[5*npmax], 1); 
    return 0;
}
