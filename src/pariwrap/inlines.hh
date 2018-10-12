inline pari p_Fl_to_Flx(ulong x, long sv)
{ return Fl_to_Flx(x, sv); }
inline pari Flc_to_ZC(pari::tmp z)
{ return Flc_to_ZC(z.g); }
inline pari Flm_to_FlxV(pari::tmp x, long sv)
{ return Flm_to_FlxV(x.g, sv); }
inline pari Flm_to_FlxX(pari::tmp x, long v, long w)
{ return Flm_to_FlxX(x.g, v, w); }
inline pari Flm_to_ZM(pari::tmp z)
{ return Flm_to_ZM(z.g); }
inline pari Flv_to_Flx(pari::tmp x, long vs)
{ return Flv_to_Flx(x.g, vs); }
inline pari Flv_to_ZV(pari::tmp z)
{ return Flv_to_ZV(z.g); }
inline pari Flv_polint(pari::tmp xa, pari::tmp ya, ulong p, long vs)
{ return Flv_polint(xa.g, ya.g, p, vs); }
inline pari Flv_roots_to_pol(pari::tmp a, ulong p, long vs)
{ return Flv_roots_to_pol(a.g, p, vs); }
inline pari Flx_Fl_mul(pari::tmp y, ulong x, ulong p)
{ return Flx_Fl_mul(y.g, x, p); }
inline pari Flx_to_Flv(pari::tmp x, long N)
{ return Flx_to_Flv(x.g, N); }
inline pari Flx_to_ZX(pari::tmp z)
{ return Flx_to_ZX(z.g); }
inline pari Flx_to_ZX_inplace(pari::tmp z)
{ return Flx_to_ZX_inplace(z.g); }
inline pari Flx_add(pari::tmp x, pari::tmp y, ulong p)
{ return Flx_add(x.g, y.g, p); }
inline pari Flx_deriv(pari::tmp z, ulong p)
{ return Flx_deriv(z.g, p); }
inline pari Flx_div_by_X_x(pari::tmp a, ulong x, ulong p, ulong * rem)
{ return Flx_div_by_X_x(a.g, x, p, rem); }
inline ulong Flx_eval(pari::tmp x, ulong y, ulong p)
{ return Flx_eval(x.g, y, p); }
inline pari Flx_gcd(pari::tmp a, pari::tmp b, ulong p)
{ return Flx_gcd(a.g, b.g, p); }
inline pari Flx_gcd_i(pari::tmp a, pari::tmp b, ulong p)
{ return Flx_gcd_i(a.g, b.g, p); }
inline pari Flx_invmontgomery(pari::tmp T, ulong p)
{ return Flx_invmontgomery(T.g, p); }
inline int Flx_is_squarefree(pari::tmp z, ulong p)
{ return Flx_is_squarefree(z.g, p); }
inline pari Flx_mul(pari::tmp x, pari::tmp y, ulong p)
{ return Flx_mul(x.g, y.g, p); }
inline pari Flx_neg(pari::tmp x, ulong p)
{ return Flx_neg(x.g, p); }
inline pari Flx_neg_inplace(pari::tmp x, ulong p)
{ return Flx_neg_inplace(x.g, p); }
inline pari Flx_normalize(pari::tmp z, ulong p)
{ return Flx_normalize(z.g, p); }
inline pari Flx_pow(pari::tmp x, long n, ulong p)
{ return Flx_pow(x.g, n, p); }
inline pari Flx_recip(pari::tmp x)
{ return Flx_recip(x.g); }
inline pari Flx_red(pari::tmp z, ulong p)
{ return Flx_red(z.g, p); }
inline pari Flx_rem_montgomery(pari::tmp x, pari::tmp mg, pari::tmp T, ulong p)
{ return Flx_rem_montgomery(x.g, mg.g, T.g, p); }
inline pari Flx_rem(pari::tmp x, pari::tmp y, ulong p)
{ return Flx_rem(x.g, y.g, p); }
inline pari Flx_renormalize(pari::tmp x, long l)
{ return Flx_renormalize(x.g, l); }
inline ulong Flx_resultant(pari::tmp a, pari::tmp b, ulong p)
{ return Flx_resultant(a.g, b.g, p); }
inline pari Flx_shift(pari::tmp a, long n)
{ return Flx_shift(a.g, n); }
inline pari Flx_sqr(pari::tmp x, ulong p)
{ return Flx_sqr(x.g, p); }
inline pari Flx_sub(pari::tmp x, pari::tmp y, ulong p)
{ return Flx_sub(x.g, y.g, p); }
inline long Flx_valuation(pari::tmp x)
{ return Flx_valuation(x.g); }
inline pari FlxC_to_ZXC(pari::tmp x)
{ return FlxC_to_ZXC(x.g); }
inline pari FlxM_to_ZXM(pari::tmp z)
{ return FlxM_to_ZXM(z.g); }
inline pari FlxV_Flc_mul(pari::tmp V, pari::tmp W, ulong p)
{ return FlxV_Flc_mul(V.g, W.g, p); }
inline pari FlxV_to_Flm(pari::tmp v, long n)
{ return FlxV_to_Flm(v.g, n); }
inline pari FlxX_add(pari::tmp P, pari::tmp Q, ulong p)
{ return FlxX_add(P.g, Q.g, p); }
inline pari FlxX_renormalize(pari::tmp x, long lx)
{ return FlxX_renormalize(x.g, lx); }
inline pari FlxX_shift(pari::tmp a, long n)
{ return FlxX_shift(a.g, n); }
inline pari FlxX_to_Flm(pari::tmp v, long n)
{ return FlxX_to_Flm(v.g, n); }
inline pari FlxX_to_ZXX(pari::tmp B)
{ return FlxX_to_ZXX(B.g); }
inline pari FlxYqQ_pow(pari::tmp x, pari::tmp n, pari::tmp S, pari::tmp T, ulong p)
{ return FlxYqQ_pow(x.g, n.g, S.g, T.g, p); }
inline pari Flxq_inv(pari::tmp x, pari::tmp T, ulong p)
{ return Flxq_inv(x.g, T.g, p); }
inline pari Flxq_invsafe(pari::tmp x, pari::tmp T, ulong p)
{ return Flxq_invsafe(x.g, T.g, p); }
inline pari Flxq_mul(pari::tmp y, pari::tmp x, pari::tmp T, ulong p)
{ return Flxq_mul(y.g, x.g, T.g, p); }
inline pari Flxq_pow(pari::tmp x, pari::tmp n, pari::tmp T, ulong p)
{ return Flxq_pow(x.g, n.g, T.g, p); }
inline pari Flxq_powers(pari::tmp x, long l, pari::tmp T, ulong p)
{ return Flxq_powers(x.g, l, T.g, p); }
inline pari Flxq_sqr(pari::tmp y, pari::tmp T, ulong p)
{ return Flxq_sqr(y.g, T.g, p); }
inline pari FlxqV_roots_to_pol(pari::tmp V, pari::tmp T, ulong p, long v)
{ return FlxqV_roots_to_pol(V.g, T.g, p, v); }
inline pari FlxqX_normalize(pari::tmp z, pari::tmp T, ulong p)
{ return FlxqX_normalize(z.g, T.g, p); }
inline pari FlxqX_Flxq_mul(pari::tmp P, pari::tmp U, pari::tmp T, ulong p)
{ return FlxqX_Flxq_mul(P.g, U.g, T.g, p); }
inline pari FlxqX_red(pari::tmp z, pari::tmp T, ulong p)
{ return FlxqX_red(z.g, T.g, p); }
inline pari FlxqX_mul(pari::tmp x, pari::tmp y, pari::tmp T, ulong p)
{ return FlxqX_mul(x.g, y.g, T.g, p); }
inline pari FlxqX_safegcd(pari::tmp P, pari::tmp Q, pari::tmp T, ulong p)
{ return FlxqX_safegcd(P.g, Q.g, T.g, p); }
inline pari FlxqX_sqr(pari::tmp x, pari::tmp T, ulong p)
{ return FlxqX_sqr(x.g, T.g, p); }
inline pari FlxqXQ_pow(pari::tmp x, pari::tmp n, pari::tmp S, pari::tmp T, ulong p)
{ return FlxqXQ_pow(x.g, n.g, S.g, T.g, p); }
inline pari FlxqXV_prod(pari::tmp V, pari::tmp T, ulong p)
{ return FlxqXV_prod(V.g, T.g, p); }
inline pari Z_to_Flx(pari::tmp x, ulong p, long v)
{ return Z_to_Flx(x.g, p, v); }
inline pari ZM_to_Flm(pari::tmp x, ulong p)
{ return ZM_to_Flm(x.g, p); }
inline pari ZV_to_Flv(pari::tmp x, ulong p)
{ return ZV_to_Flv(x.g, p); }
inline pari ZX_to_Flx(pari::tmp x, ulong p)
{ return ZX_to_Flx(x.g, p); }
inline pari ZXV_to_FlxV(pari::tmp v, ulong p)
{ return ZXV_to_FlxV(v.g, p); }
inline pari ZXX_to_FlxX(pari::tmp B, ulong p, long v)
{ return ZXX_to_FlxX(B.g, p, v); }
inline pari ZXXV_to_FlxXV(pari::tmp V, ulong p, long v)
{ return ZXXV_to_FlxXV(V.g, p, v); }
inline pari p_polx_Flx(long sv)
{ return polx_Flx(sv); }
inline pari p_zero_Flx(long sv)
{ return zero_Flx(sv); }
inline pari Flm_Flc_mul(pari::tmp x, pari::tmp y, ulong p)
{ return Flm_Flc_mul(x.g, y.g, p); }
inline pari Flm_deplin(pari::tmp x, ulong p)
{ return Flm_deplin(x.g, p); }
inline pari Flm_gauss(pari::tmp a, pari::tmp b, ulong p)
{ return Flm_gauss(a.g, b.g, p); }
inline pari Flm_indexrank(pari::tmp x, ulong p)
{ return Flm_indexrank(x.g, p); }
inline pari Flm_inv(pari::tmp x, ulong p)
{ return Flm_inv(x.g, p); }
inline pari Flm_ker(pari::tmp x, ulong p)
{ return Flm_ker(x.g, p); }
inline pari Flm_ker_sp(pari::tmp x, ulong p, long deplin)
{ return Flm_ker_sp(x.g, p, deplin); }
inline pari Flm_mul(pari::tmp x, pari::tmp y, ulong p)
{ return Flm_mul(x.g, y.g, p); }
inline pari FlxqM_ker(pari::tmp x, pari::tmp T, ulong p)
{ return FlxqM_ker(x.g, T.g, p); }
inline pari FpC_Fp_mul(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpC_Fp_mul(x.g, y.g, p.g); }
inline pari FpC_FpV_mul(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpC_FpV_mul(x.g, y.g, p.g); }
inline pari FpM_FpC_mul(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpM_FpC_mul(x.g, y.g, p.g); }
inline pari FpM_deplin(pari::tmp x, pari::tmp p)
{ return FpM_deplin(x.g, p.g); }
inline pari FpM_gauss(pari::tmp a, pari::tmp b, pari::tmp p)
{ return FpM_gauss(a.g, b.g, p.g); }
inline pari FpM_image(pari::tmp x, pari::tmp p)
{ return FpM_image(x.g, p.g); }
inline pari FpM_intersect(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpM_intersect(x.g, y.g, p.g); }
inline pari FpM_inv(pari::tmp x, pari::tmp p)
{ return FpM_inv(x.g, p.g); }
inline pari FpM_invimage(pari::tmp m, pari::tmp v, pari::tmp p)
{ return FpM_invimage(m.g, v.g, p.g); }
inline pari FpM_ker(pari::tmp x, pari::tmp p)
{ return FpM_ker(x.g, p.g); }
inline pari FpM_mul(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpM_mul(x.g, y.g, p.g); }
inline long FpM_rank(pari::tmp x, pari::tmp p)
{ return FpM_rank(x.g, p.g); }
inline pari FpM_indexrank(pari::tmp x, pari::tmp p)
{ return FpM_indexrank(x.g, p.g); }
inline pari FpM_suppl(pari::tmp x, pari::tmp p)
{ return FpM_suppl(x.g, p.g); }
inline pari FpV_FpC_mul(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpV_FpC_mul(x.g, y.g, p.g); }
inline pari FqM_gauss(pari::tmp a, pari::tmp b, pari::tmp T, pari::tmp p)
{ return FqM_gauss(a.g, b.g, T.g, p.g); }
inline pari FqM_ker(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FqM_ker(x.g, T.g, p.g); }
inline pari FqM_suppl(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FqM_suppl(x.g, T.g, p.g); }
inline pari QM_inv(pari::tmp M, pari::tmp dM)
{ return QM_inv(M.g, dM.g); }
inline pari ZM_inv(pari::tmp M, pari::tmp dM)
{ return ZM_inv(M.g, dM.g); }
inline void appendL(pari::tmp x, pari::tmp t)
{ appendL(x.g, t.g); }
inline pari p_cget1(long l, long t)
{ return cget1(l, t); }
inline pari concat(pari::tmp x, pari::tmp y)
{ return concat(x.g, y.g); }
inline pari shallowconcat(pari::tmp x, pari::tmp y)
{ return shallowconcat(x.g, y.g); }
inline pari concatsp3(pari::tmp x, pari::tmp y, pari::tmp z)
{ return concatsp3(x.g, y.g, z.g); }
inline pari deplin(pari::tmp x)
{ return deplin(x.g); }
inline pari det(pari::tmp a)
{ return det(a.g); }
inline pari det0(pari::tmp a, long flag)
{ return det0(a.g, flag); }
inline pari det2(pari::tmp a)
{ return det2(a.g); }
inline pari detint(pari::tmp x)
{ return detint(x.g); }
inline pari diagonal(pari::tmp x)
{ return diagonal(x.g); }
inline pari eigen(pari::tmp x, long prc = prec)
{ return eigen(x.g, prc); }
inline pari extract(pari::tmp x, pari::tmp l)
{ return extract(x.g, l.g); }
inline pari extract0(pari::tmp x, pari::tmp l1, pari::tmp l2)
{ return extract0(x.g, l1.g, l2.g); }
inline pari gaddmat(pari::tmp x, pari::tmp y)
{ return gaddmat(x.g, y.g); }
inline pari gaddmat_i(pari::tmp x, pari::tmp y)
{ return gaddmat_i(x.g, y.g); }
inline pari gauss(pari::tmp a, pari::tmp b)
{ return gauss(a.g, b.g); }
inline pari gaussmodulo(pari::tmp M, pari::tmp D, pari::tmp Y)
{ return gaussmodulo(M.g, D.g, Y.g); }
inline pari gaussmodulo2(pari::tmp M, pari::tmp D, pari::tmp Y)
{ return gaussmodulo2(M.g, D.g, Y.g); }
inline pari gscalcol(pari::tmp x, long n)
{ return gscalcol(x.g, n); }
inline pari gscalcol_i(pari::tmp x, long n)
{ return gscalcol_i(x.g, n); }
inline pari gscalmat(pari::tmp x, long n)
{ return gscalmat(x.g, n); }
inline pari p_gscalsmat(long x, long n)
{ return gscalsmat(x, n); }
inline pari gtomat(pari::tmp x)
{ return gtomat(x.g); }
inline pari gtrans(pari::tmp x)
{ return gtrans(x.g); }
inline pari shallowtrans(pari::tmp x)
{ return shallowtrans(x.g); }
inline int hnfdivide(pari::tmp A, pari::tmp B)
{ return hnfdivide(A.g, B.g); }
inline pari p_matid(long n)
{ return matid(n); }
inline pari p_matid_Flm(long n)
{ return matid_Flm(n); }
inline pari image(pari::tmp x)
{ return image(x.g); }
inline pari image2(pari::tmp x)
{ return image2(x.g); }
inline pari imagecompl(pari::tmp x)
{ return imagecompl(x.g); }
inline pari indexrank(pari::tmp x)
{ return indexrank(x.g); }
inline pari inverseimage(pari::tmp mat, pari::tmp y)
{ return inverseimage(mat.g, y.g); }
inline long isdiagonal(pari::tmp x)
{ return isdiagonal(x.g); }
inline long isscalarmat(pari::tmp x, pari::tmp s)
{ return isscalarmat(x.g, s.g); }
inline pari ker(pari::tmp x)
{ return ker(x.g); }
inline pari keri(pari::tmp x)
{ return keri(x.g); }
inline pari matextract(pari::tmp x, pari::tmp l1, pari::tmp l2)
{ return matextract(x.g, l1.g, l2.g); }
inline pari matimage0(pari::tmp x, long flag)
{ return matimage0(x.g, flag); }
inline pari matker0(pari::tmp x, long flag)
{ return matker0(x.g, flag); }
inline pari matmuldiagonal(pari::tmp x, pari::tmp d)
{ return matmuldiagonal(x.g, d.g); }
inline pari matmultodiagonal(pari::tmp x, pari::tmp y)
{ return matmultodiagonal(x.g, y.g); }
inline pari matsolvemod0(pari::tmp M, pari::tmp D, pari::tmp Y, long flag)
{ return matsolvemod0(M.g, D.g, Y.g, flag); }
inline pari mattodiagonal(pari::tmp m)
{ return mattodiagonal(m.g); }
inline pari mattodiagonal_i(pari::tmp m)
{ return mattodiagonal_i(m.g); }
//inline long rank(pari::tmp x)
//{ return rank(x.g); }
inline pari row(pari::tmp A, long x1)
{ return row(A.g, x1); }
inline pari row_i(pari::tmp A, long x0, long x1, long x2)
{ return row_i(A.g, x0, x1, x2); }
inline pari rowcopy(pari::tmp A, long x0)
{ return rowcopy(A.g, x0); }
inline pari rowslice(pari::tmp A, long x1, long x2)
{ return rowslice(A.g, x1, x2); }
inline pari rowslicepermute(pari::tmp A, pari::tmp p, long x1, long x2)
{ return rowslicepermute(A.g, p.g, x1, x2); }
inline pari rowpermute(pari::tmp A, pari::tmp p)
{ return rowpermute(A.g, p.g); }
inline pari sindexrank(pari::tmp x)
{ return sindexrank(x.g); }
inline pari sum(pari::tmp v, long a, long b)
{ return sum(v.g, a, b); }
inline pari suppl(pari::tmp x)
{ return suppl(x.g); }
inline pari vconcat(pari::tmp A, pari::tmp B)
{ return vconcat(A.g, B.g); }
inline pari vecslice(pari::tmp A, long y1, long y2)
{ return vecslice(A.g, y1, y2); }
inline pari vecslicepermute(pari::tmp A, pari::tmp p, long y1, long y2)
{ return vecslicepermute(A.g, p.g, y1, y2); }
inline pari vecpermute(pari::tmp A, pari::tmp p)
{ return vecpermute(A.g, p.g); }
inline pari QuickNormL1(pari::tmp x, long prc = prec)
{ return QuickNormL1(x.g, prc); }
inline pari QuickNormL2(pari::tmp x, long prc = prec)
{ return QuickNormL2(x.g, prc); }
inline int RgM_ishnf(pari::tmp x)
{ return RgM_ishnf(x.g); }
inline pari RgXQ_norm(pari::tmp x, pari::tmp T)
{ return RgXQ_norm(x.g, T.g); }
inline pari ZM_to_zm(pari::tmp z)
{ return ZM_to_zm(z.g); }
inline int ZM_ishnf(pari::tmp x)
{ return ZM_ishnf(x.g); }
inline pari ZV_add(pari::tmp x, pari::tmp y)
{ return ZV_add(x.g, y.g); }
inline pari ZV_sub(pari::tmp x, pari::tmp y)
{ return ZV_sub(x.g, y.g); }
inline pari ZV_to_nv(pari::tmp z)
{ return ZV_to_nv(z.g); }
inline pari adj(pari::tmp x)
{ return adj(x.g); }
inline pari assmat(pari::tmp x)
{ return assmat(x.g); }
inline pari caract(pari::tmp x, long v)
{ return caract(x.g, v); }
inline pari caract2(pari::tmp p, pari::tmp x, long v)
{ return caract2(p.g, x.g, v); }
inline pari caradj0(pari::tmp x, long v)
{ return caradj0(x.g, v); }
inline pari carhess(pari::tmp x, long v)
{ return carhess(x.g, v); }
inline pari charpoly0(pari::tmp x, long v, long flag)
{ return charpoly0(x.g, v, flag); }
inline pari conjvec(pari::tmp x, long prc = prec)
{ return conjvec(x.g, prc); }
inline pari gconj(pari::tmp x)
{ return gconj(x.g); }
inline pari gnorm(pari::tmp x)
{ return gnorm(x.g); }
inline pari gnorml1(pari::tmp x, long prc = prec)
{ return gnorml1(x.g, prc); }
inline pari gnorml2(pari::tmp x)
{ return gnorml2(x.g); }
inline pari gsmith(pari::tmp x)
{ return gsmith(x.g); }
inline pari gsmith2(pari::tmp x)
{ return gsmith2(x.g); }
inline pari gtrace(pari::tmp x)
{ return gtrace(x.g); }
inline pari hess(pari::tmp x)
{ return hess(x.g); }
inline pari hnf(pari::tmp x)
{ return hnf(x.g); }
inline pari hnfall(pari::tmp x)
{ return hnfall(x.g); }
inline pari hnflll(pari::tmp x)
{ return hnflll(x.g); }
inline pari hnfmod(pari::tmp x, pari::tmp detmat)
{ return hnfmod(x.g, detmat.g); }
inline pari hnfmodid(pari::tmp x, pari::tmp p)
{ return hnfmodid(x.g, p.g); }
inline pari hnfmodidpart(pari::tmp x, pari::tmp p)
{ return hnfmodidpart(x.g, p.g); }
inline pari hnfperm(pari::tmp x)
{ return hnfperm(x.g); }
inline pari intersect(pari::tmp x, pari::tmp y)
{ return intersect(x.g, y.g); }
inline pari jacobi(pari::tmp a, long prc = prec)
{ return jacobi(a.g, prc); }
inline pari matfrobenius(pari::tmp M, long flag, long v)
{ return matfrobenius(M.g, flag, v); }
inline pari mathnf0(pari::tmp x, long flag)
{ return mathnf0(x.g, flag); }
inline pari matrixqz(pari::tmp x, pari::tmp pp)
{ return matrixqz(x.g, pp.g); }
inline pari matrixqz0(pari::tmp x, pari::tmp pp)
{ return matrixqz0(x.g, pp.g); }
inline pari matrixqz2(pari::tmp x)
{ return matrixqz2(x.g); }
inline pari matrixqz3(pari::tmp x)
{ return matrixqz3(x.g); }
inline pari matsnf0(pari::tmp x, long flag)
{ return matsnf0(x.g, flag); }
inline pari minpoly(pari::tmp x, long v)
{ return minpoly(x.g, v); }
inline pari signat(pari::tmp a)
{ return signat(a.g); }
inline pari smith(pari::tmp x)
{ return smith(x.g); }
inline pari smith2(pari::tmp x)
{ return smith2(x.g); }
inline pari smithclean(pari::tmp z)
{ return smithclean(z.g); }
inline pari sqred(pari::tmp a)
{ return sqred(a.g); }
inline pari sqred1(pari::tmp a)
{ return sqred1(a.g); }
inline pari sqred3(pari::tmp a)
{ return sqred3(a.g); }
inline pari zm_to_ZM(pari::tmp z)
{ return zm_to_ZM(z.g); }
inline pari zx_to_ZX(pari::tmp z)
{ return zx_to_ZX(z.g); }
inline pari p_gp_read_str(char * t)
{ return gp_read_str(t); }
inline pari p_readseq(char * t)
{ return readseq(t); }
inline pari p_strtoGENstr(const char * s)
{ return strtoGENstr(s); }
inline pari p_strtoi(char * s)
{ return strtoi(s); }
inline pari p_strtor(char * s, long prc = prec)
{ return strtor(s, prc); }
inline pari type0(pari::tmp x)
{ return type0(x.g); }
inline long isprimeAPRCL(pari::tmp N)
{ return isprimeAPRCL(N.g); }
inline void check_quaddisc(pari::tmp x, long * s, long * r, char * f)
{ check_quaddisc(x.g, s, r, f); }
inline void check_quaddisc_real(pari::tmp x, long * r, char * f)
{ check_quaddisc_real(x.g, r, f); }
inline void check_quaddisc_imag(pari::tmp x, long * r, char * f)
{ check_quaddisc_imag(x.g, r, f); }
inline pari compimag(pari::tmp x, pari::tmp y)
{ return compimag(x.g, y.g); }
inline pari compimagraw(pari::tmp x, pari::tmp y)
{ return compimagraw(x.g, y.g); }
inline pari compraw(pari::tmp x, pari::tmp y)
{ return compraw(x.g, y.g); }
inline pari compreal(pari::tmp x, pari::tmp y)
{ return compreal(x.g, y.g); }
inline pari comprealraw(pari::tmp x, pari::tmp y)
{ return comprealraw(x.g, y.g); }
inline pari nucomp(pari::tmp x, pari::tmp y, pari::tmp l)
{ return nucomp(x.g, y.g, l.g); }
inline pari nudupl(pari::tmp x, pari::tmp l)
{ return nudupl(x.g, l.g); }
inline pari nupow(pari::tmp x, pari::tmp n)
{ return nupow(x.g, n.g); }
inline pari powraw(pari::tmp x, long n)
{ return powraw(x.g, n); }
inline pari powrealraw(pari::tmp x, long n)
{ return powrealraw(x.g, n); }
inline pari primeform(pari::tmp x, pari::tmp p, long prc = prec)
{ return primeform(x.g, p.g, prc); }
inline pari Qfb0(pari::tmp x, pari::tmp y, pari::tmp z, pari::tmp d, long prc = prec)
{ return Qfb0(x.g, y.g, z.g, d.g, prc); }
inline pari qfbimagsolvep(pari::tmp Q, pari::tmp p)
{ return qfbimagsolvep(Q.g, p.g); }
inline pari qfbrealsolvep(pari::tmp Q, pari::tmp p)
{ return qfbrealsolvep(Q.g, p.g); }
inline pari qfbred0(pari::tmp x, long flag, pari::tmp D, pari::tmp isqrtD, pari::tmp sqrtD)
{ return qfbred0(x.g, flag, D.g, isqrtD.g, sqrtD.g); }
inline pari qfbsolve(pari::tmp Q, pari::tmp n)
{ return qfbsolve(Q.g, n.g); }
inline pari qfi(pari::tmp x, pari::tmp y, pari::tmp z)
{ return qfi(x.g, y.g, z.g); }
inline pari qfr(pari::tmp x, pari::tmp y, pari::tmp z, pari::tmp d)
{ return qfr(x.g, y.g, z.g, d.g); }
inline pari quadgen(pari::tmp x)
{ return quadgen(x.g); }
inline pari quadpoly(pari::tmp x)
{ return quadpoly(x.g); }
inline pari quadpoly0(pari::tmp x, long v)
{ return quadpoly0(x.g, v); }
inline pari redimag(pari::tmp x)
{ return redimag(x.g); }
inline pari redreal(pari::tmp x)
{ return redreal(x.g); }
inline pari redrealnod(pari::tmp x, pari::tmp isqrtD)
{ return redrealnod(x.g, isqrtD.g); }
inline pari rhoreal(pari::tmp x)
{ return rhoreal(x.g); }
inline pari rhorealnod(pari::tmp x, pari::tmp isqrtD)
{ return rhorealnod(x.g, isqrtD.g); }
inline pari sqcompimag(pari::tmp x)
{ return sqcompimag(x.g); }
inline pari sqcompreal(pari::tmp x)
{ return sqcompreal(x.g); }
inline pari Fp_inv(pari::tmp a, pari::tmp m)
{ return Fp_inv(a.g, m.g); }
inline pari Fp_invsafe(pari::tmp a, pari::tmp m)
{ return Fp_invsafe(a.g, m.g); }
inline pari Fp_pow(pari::tmp a, pari::tmp n, pari::tmp m)
{ return Fp_pow(a.g, n.g, m.g); }
inline pari Fp_pows(pari::tmp A, long k, pari::tmp N)
{ return Fp_pows(A.g, k, N.g); }
inline pari Fp_powu(pari::tmp x, ulong k, pari::tmp p)
{ return Fp_powu(x.g, k, p.g); }
inline pari Fp_sqrt(pari::tmp a, pari::tmp p)
{ return Fp_sqrt(a.g, p.g); }
inline pari bestappr0(pari::tmp x, pari::tmp a, pari::tmp b)
{ return bestappr0(x.g, a.g, b.g); }
inline pari bestappr(pari::tmp x, pari::tmp k)
{ return bestappr(x.g, k.g); }
inline pari chinese1(pari::tmp x)
{ return chinese1(x.g); }
inline pari chinese(pari::tmp x, pari::tmp y)
{ return chinese(x.g, y.g); }
inline pari classno2(pari::tmp x)
{ return classno2(x.g); }
inline pari classno(pari::tmp x)
{ return classno(x.g); }
inline pari contfrac0(pari::tmp x, pari::tmp b, long flag)
{ return contfrac0(x.g, b.g, flag); }
inline pari p_fibo(long n)
{ return fibo(n); }
inline pari fundunit(pari::tmp x)
{ return fundunit(x.g); }
inline pari gboundcf(pari::tmp x, long k)
{ return gboundcf(x.g, k); }
inline pari gcf2(pari::tmp b, pari::tmp x)
{ return gcf2(b.g, x.g); }
inline pari gcf(pari::tmp x)
{ return gcf(x.g); }
inline pari gener(pari::tmp m)
{ return gener(m.g); }
inline ulong gener_Fl_local(ulong p, pari::tmp L)
{ return gener_Fl_local(p, L.g); }
inline pari gener_Fp_local(pari::tmp p, pari::tmp L)
{ return gener_Fp_local(p.g, L.g); }
inline pari gener_Fp(pari::tmp p)
{ return gener_Fp(p.g); }
inline pari gfundunit(pari::tmp x)
{ return gfundunit(x.g); }
inline pari ggener(pari::tmp m)
{ return ggener(m.g); }
inline pari gisfundamental(pari::tmp x)
{ return gisfundamental(x.g); }
inline pari gisprime(pari::tmp x, long flag)
{ return gisprime(x.g, flag); }
inline pari gispseudoprime(pari::tmp x, long flag)
{ return gispseudoprime(x.g, flag); }
inline pari gispsp(pari::tmp x)
{ return gispsp(x.g); }
inline pari gissquare(pari::tmp x)
{ return gissquare(x.g); }
inline pari gkrogs(pari::tmp x, long y)
{ return gkrogs(x.g, y); }
inline pari gkronecker(pari::tmp x, pari::tmp y)
{ return gkronecker(x.g, y.g); }
inline pari gmillerrabin(pari::tmp n, long k)
{ return gmillerrabin(n.g, k); }
inline pari gnextprime(pari::tmp n)
{ return gnextprime(n.g); }
inline pari gprecprime(pari::tmp n)
{ return gprecprime(n.g); }
inline pari gracine(pari::tmp a)
{ return gracine(a.g); }
inline pari gregula(pari::tmp x, long prc = prec)
{ return gregula(x.g, prc); }
inline pari hclassno(pari::tmp x)
{ return hclassno(x.g); }
inline long hil0(pari::tmp x, pari::tmp y, pari::tmp p)
{ return hil0(x.g, y.g, p.g); }
inline long hil(pari::tmp x, pari::tmp y, pari::tmp p)
{ return hil(x.g, y.g, p.g); }
inline long isfundamental(pari::tmp x)
{ return isfundamental(x.g); }
inline long isprime(pari::tmp x)
{ return isprime(x.g); }
inline long isprimeSelfridge(pari::tmp x)
{ return isprimeSelfridge(x.g); }
inline long ispseudoprime(pari::tmp x, long flag)
{ return ispseudoprime(x.g, flag); }
inline long ispsp(pari::tmp x)
{ return ispsp(x.g); }
inline long krois(pari::tmp x, long y)
{ return krois(x.g, y); }
inline long kronecker(pari::tmp x, pari::tmp y)
{ return kronecker(x.g, y.g); }
inline long krosi(long s, pari::tmp x)
{ return krosi(s, x.g); }
inline pari lcmii(pari::tmp a, pari::tmp b)
{ return lcmii(a.g, b.g); }
inline pari p_mpfact(long n)
{ return mpfact(n); }
inline pari order(pari::tmp x)
{ return order(x.g); }
inline pari pnqn(pari::tmp x)
{ return pnqn(x.g); }
inline pari qfbclassno0(pari::tmp x, long flag)
{ return qfbclassno0(x.g, flag); }
inline pari quaddisc(pari::tmp x)
{ return quaddisc(x.g); }
inline pari racine(pari::tmp a)
{ return racine(a.g); }
inline pari regula(pari::tmp x, long prc = prec)
{ return regula(x.g, prc); }
inline pari p_seq_umul(ulong a, ulong b)
{ return seq_umul(a, b); }
inline pari znorder(pari::tmp x, pari::tmp o)
{ return znorder(x.g, o.g); }
inline pari znstar(pari::tmp x)
{ return znstar(x.g); }
inline long Z_issquarefree(pari::tmp x)
{ return Z_issquarefree(x.g); }
inline pari addprimes(pari::tmp primes)
{ return addprimes(primes.g); }
inline pari auxdecomp(pari::tmp n, long all)
{ return auxdecomp(n.g, all); }
inline long bigomega(pari::tmp n)
{ return bigomega(n.g); }
inline pari binaire(pari::tmp x)
{ return binaire(x.g); }
inline long bittest(pari::tmp x, long n)
{ return bittest(x.g, n); }
inline pari boundfact(pari::tmp n, long lim)
{ return boundfact(n.g, lim); }
inline pari core(pari::tmp n)
{ return core(n.g); }
inline pari corepartial(pari::tmp n, long l)
{ return corepartial(n.g, l); }
inline pari core0(pari::tmp n, long flag)
{ return core0(n.g, flag); }
inline pari core2(pari::tmp n)
{ return core2(n.g); }
inline pari core2partial(pari::tmp n, long l)
{ return core2partial(n.g, l); }
inline pari coredisc(pari::tmp n)
{ return coredisc(n.g); }
inline pari coredisc0(pari::tmp n, long flag)
{ return coredisc0(n.g, flag); }
inline pari coredisc2(pari::tmp n)
{ return coredisc2(n.g); }
inline pari divisors(pari::tmp n)
{ return divisors(n.g); }
inline pari factorint(pari::tmp n, long flag)
{ return factorint(n.g, flag); }
inline pari p_factoru(ulong n)
{ return factoru(n); }
inline pari p_factoru_pow(ulong n)
{ return factoru_pow(n); }
inline pari gbigomega(pari::tmp n)
{ return gbigomega(n.g); }
inline pari gbitand(pari::tmp x, pari::tmp y)
{ return gbitand(x.g, y.g); }
inline pari gbitneg(pari::tmp x, long n)
{ return gbitneg(x.g, n); }
inline pari gbitnegimply(pari::tmp x, pari::tmp y)
{ return gbitnegimply(x.g, y.g); }
inline pari gbitor(pari::tmp x, pari::tmp y)
{ return gbitor(x.g, y.g); }
inline pari gbittest(pari::tmp x, pari::tmp n)
{ return gbittest(x.g, n.g); }
inline pari gbitxor(pari::tmp x, pari::tmp y)
{ return gbitxor(x.g, y.g); }
inline pari gboundfact(pari::tmp n, long lim)
{ return gboundfact(n.g, lim); }
inline pari gissquarefree(pari::tmp x)
{ return gissquarefree(x.g); }
inline pari gmu(pari::tmp n)
{ return gmu(n.g); }
inline pari gnumbdiv(pari::tmp n)
{ return gnumbdiv(n.g); }
inline pari gomega(pari::tmp n)
{ return gomega(n.g); }
inline pari gphi(pari::tmp n)
{ return gphi(n.g); }
inline pari gsumdiv(pari::tmp n)
{ return gsumdiv(n.g); }
inline pari gsumdivk(pari::tmp n, long k)
{ return gsumdivk(n.g, k); }
inline long issquarefree(pari::tmp x)
{ return issquarefree(x.g); }
inline long mu(pari::tmp n)
{ return mu(n.g); }
inline pari numbdiv(pari::tmp n)
{ return numbdiv(n.g); }
inline long omega(pari::tmp n)
{ return omega(n.g); }
inline pari phi(pari::tmp n)
{ return phi(n.g); }
inline pari p_prime(long n)
{ return prime(n); }
inline pari primepi(pari::tmp x)
{ return primepi(x.g); }
inline pari p_primes(long n)
{ return primes(n); }
inline pari removeprimes(pari::tmp primes)
{ return removeprimes(primes.g); }
inline pari smallfact(pari::tmp n)
{ return smallfact(n.g); }
inline pari sumdiv(pari::tmp n)
{ return sumdiv(n.g); }
inline pari sumdivk(pari::tmp n, long k)
{ return sumdivk(n.g, k); }
inline pari Z_factor(pari::tmp n)
{ return Z_factor(n.g); }
inline pari T2_from_embed(pari::tmp x, long r1)
{ return T2_from_embed(x.g, r1); }
inline void check_ZX(pari::tmp x, char * s)
{ check_ZX(x.g, s); }
inline void check_ZXY(pari::tmp x, char * s)
{ check_ZXY(x.g, s); }
inline pari check_units(pari::tmp x, char * f)
{ return check_units(x.g, f); }
inline void checkbid(pari::tmp bid)
{ checkbid(bid.g); }
inline pari checkbnf(pari::tmp bnf)
{ return checkbnf(bnf.g); }
inline void checkbnr(pari::tmp bnr)
{ checkbnr(bnr.g); }
inline void checkbnrgen(pari::tmp bnr)
{ checkbnrgen(bnr.g); }
inline void checkid(pari::tmp x, long N)
{ checkid(x.g, N); }
inline pari checknf(pari::tmp nf)
{ return checknf(nf.g); }
inline pari checknfelt_mod(pari::tmp nf, pari::tmp x, char * s)
{ return checknfelt_mod(nf.g, x.g, s); }
inline void checkprimeid(pari::tmp bid)
{ checkprimeid(bid.g); }
inline void checkrnf(pari::tmp rnf)
{ checkrnf(rnf.g); }
inline pari dirzetak(pari::tmp nf, pari::tmp b)
{ return dirzetak(nf.g, b.g); }
inline pari factoredpolred(pari::tmp x, pari::tmp fa)
{ return factoredpolred(x.g, fa.g); }
inline pari factoredpolred2(pari::tmp x, pari::tmp fa)
{ return factoredpolred2(x.g, fa.g); }
inline pari galois(pari::tmp x, long prc = prec)
{ return galois(x.g, prc); }
inline pari galoisapply(pari::tmp nf, pari::tmp aut, pari::tmp x)
{ return galoisapply(nf.g, aut.g, x.g); }
inline pari get_bnf(pari::tmp x, long * t)
{ return get_bnf(x.g, t); }
inline pari get_nf(pari::tmp x, long * t)
{ return get_nf(x.g, t); }
inline pari get_primeid(pari::tmp x)
{ return get_primeid(x.g); }
inline pari glambdak(pari::tmp nfz, pari::tmp s, long prc = prec)
{ return glambdak(nfz.g, s.g, prc); }
inline int gpolcomp(pari::tmp p1, pari::tmp p2)
{ return gpolcomp(p1.g, p2.g); }
inline pari gzetak(pari::tmp nfz, pari::tmp s, long prc = prec)
{ return gzetak(nfz.g, s.g, prc); }
inline pari gzetakall(pari::tmp nfz, pari::tmp s, long flag, long prc = prec)
{ return gzetakall(nfz.g, s.g, flag, prc); }
inline pari initalg(pari::tmp x, long prc = prec)
{ return initalg(x.g, prc); }
inline pari initalgred(pari::tmp x, long prc = prec)
{ return initalgred(x.g, prc); }
inline pari initalgred2(pari::tmp x, long prc = prec)
{ return initalgred2(x.g, prc); }
inline pari initzeta(pari::tmp pol, long prc = prec)
{ return initzeta(pol.g, prc); }
inline long nf_get_r1(pari::tmp nf)
{ return nf_get_r1(nf.g); }
inline long nf_get_r2(pari::tmp nf)
{ return nf_get_r2(nf.g); }
inline void nf_get_sign(pari::tmp nf, long * r1, long * r2)
{ nf_get_sign(nf.g, r1, r2); }
inline long nfgetprec(pari::tmp x)
{ return nfgetprec(x.g); }
inline pari nfinit0(pari::tmp x, long flag, long prc = prec)
{ return nfinit0(x.g, flag, prc); }
inline pari nfisincl(pari::tmp a, pari::tmp b)
{ return nfisincl(a.g, b.g); }
inline pari nfisisom(pari::tmp a, pari::tmp b)
{ return nfisisom(a.g, b.g); }
inline pari nfnewprec(pari::tmp nf, long prc = prec)
{ return nfnewprec(nf.g, prc); }
inline pari nfnewprec_i(pari::tmp nf, long prc = prec)
{ return nfnewprec_i(nf.g, prc); }
inline pari ordred(pari::tmp x)
{ return ordred(x.g); }
inline pari polgalois(pari::tmp x, long prc = prec)
{ return polgalois(x.g, prc); }
inline pari polred(pari::tmp x)
{ return polred(x.g); }
inline pari polred0(pari::tmp x, long flag, pari::tmp p)
{ return polred0(x.g, flag, p.g); }
inline pari polred2(pari::tmp x)
{ return polred2(x.g); }
inline pari polredabs(pari::tmp x)
{ return polredabs(x.g); }
inline pari polredabs0(pari::tmp x, long flag)
{ return polredabs0(x.g, flag); }
inline pari polredabs2(pari::tmp x)
{ return polredabs2(x.g); }
inline pari polredabsall(pari::tmp x, long flun)
{ return polredabsall(x.g, flun); }
inline pari rootsof1(pari::tmp x)
{ return rootsof1(x.g); }
inline pari smallpolred(pari::tmp x)
{ return smallpolred(x.g); }
inline pari smallpolred2(pari::tmp x)
{ return smallpolred2(x.g); }
inline pari tschirnhaus(pari::tmp x)
{ return tschirnhaus(x.g); }
inline void checkmodpr(pari::tmp modpr)
{ checkmodpr(modpr.g); }
inline pari compositum(pari::tmp pol1, pari::tmp pol2)
{ return compositum(pol1.g, pol2.g); }
inline pari compositum2(pari::tmp pol1, pari::tmp pol2)
{ return compositum2(pol1.g, pol2.g); }
inline pari discf(pari::tmp x)
{ return discf(x.g); }
inline pari discf2(pari::tmp x)
{ return discf2(x.g); }
inline pari factoreddiscf(pari::tmp x, pari::tmp p)
{ return factoreddiscf(x.g, p.g); }
inline pari ff_to_nf(pari::tmp x, pari::tmp modpr)
{ return ff_to_nf(x.g, modpr.g); }
inline pari fix_relative_pol(pari::tmp nf, pari::tmp x, int chk_lead)
{ return fix_relative_pol(nf.g, x.g, chk_lead); }
inline pari gcdpm(pari::tmp f1, pari::tmp f2, pari::tmp pm)
{ return gcdpm(f1.g, f2.g, pm.g); }
inline pari indexpartial(pari::tmp P, pari::tmp DP)
{ return indexpartial(P.g, DP.g); }
inline pari modprX(pari::tmp x, pari::tmp nf, pari::tmp modpr)
{ return modprX(x.g, nf.g, modpr.g); }
inline pari modprX_lift(pari::tmp x, pari::tmp modpr)
{ return modprX_lift(x.g, modpr.g); }
inline pari modprM(pari::tmp z, pari::tmp nf, pari::tmp modpr)
{ return modprM(z.g, nf.g, modpr.g); }
inline pari modprM_lift(pari::tmp z, pari::tmp modpr)
{ return modprM_lift(z.g, modpr.g); }
inline pari nf_to_ff(pari::tmp nf, pari::tmp x, pari::tmp modpr)
{ return nf_to_ff(nf.g, x.g, modpr.g); }
inline pari nfbasis0(pari::tmp x, long flag, pari::tmp p)
{ return nfbasis0(x.g, flag, p.g); }
inline pari nfdiscf0(pari::tmp x, long flag, pari::tmp p)
{ return nfdiscf0(x.g, flag, p.g); }
inline pari nfmodprinit(pari::tmp nf, pari::tmp pr)
{ return nfmodprinit(nf.g, pr.g); }
inline pari nfreducemodpr(pari::tmp nf, pari::tmp x, pari::tmp modpr)
{ return nfreducemodpr(nf.g, x.g, modpr.g); }
inline pari polcompositum0(pari::tmp pol1, pari::tmp pol2, long flag)
{ return polcompositum0(pol1.g, pol2.g, flag); }
inline pari primedec(pari::tmp nf, pari::tmp p)
{ return primedec(nf.g, p.g); }
inline pari rnfbasis(pari::tmp bnf, pari::tmp order)
{ return rnfbasis(bnf.g, order.g); }
inline pari rnfdedekind(pari::tmp nf, pari::tmp T, pari::tmp pr)
{ return rnfdedekind(nf.g, T.g, pr.g); }
inline pari rnfdet(pari::tmp nf, pari::tmp order)
{ return rnfdet(nf.g, order.g); }
inline pari rnfdet2(pari::tmp nf, pari::tmp A, pari::tmp I)
{ return rnfdet2(nf.g, A.g, I.g); }
inline pari rnfdiscf(pari::tmp nf, pari::tmp pol)
{ return rnfdiscf(nf.g, pol.g); }
inline pari rnfequation(pari::tmp nf, pari::tmp pol2)
{ return rnfequation(nf.g, pol2.g); }
inline pari rnfequation0(pari::tmp nf, pari::tmp pol2, long flall)
{ return rnfequation0(nf.g, pol2.g, flall); }
inline pari rnfequation2(pari::tmp nf, pari::tmp pol)
{ return rnfequation2(nf.g, pol.g); }
inline pari rnfhnfbasis(pari::tmp bnf, pari::tmp order)
{ return rnfhnfbasis(bnf.g, order.g); }
inline long rnfisfree(pari::tmp bnf, pari::tmp order)
{ return rnfisfree(bnf.g, order.g); }
inline pari rnflllgram(pari::tmp nf, pari::tmp pol, pari::tmp order, long prc = prec)
{ return rnflllgram(nf.g, pol.g, order.g, prc); }
inline pari rnfpolred(pari::tmp nf, pari::tmp pol, long prc = prec)
{ return rnfpolred(nf.g, pol.g, prc); }
inline pari rnfpolredabs(pari::tmp nf, pari::tmp pol, long flag)
{ return rnfpolredabs(nf.g, pol.g, flag); }
inline pari rnfpseudobasis(pari::tmp nf, pari::tmp pol)
{ return rnfpseudobasis(nf.g, pol.g); }
inline pari rnfsimplifybasis(pari::tmp bnf, pari::tmp order)
{ return rnfsimplifybasis(bnf.g, order.g); }
inline pari rnfsteinitz(pari::tmp nf, pari::tmp order)
{ return rnfsteinitz(nf.g, order.g); }
inline pari smalldiscf(pari::tmp x)
{ return smalldiscf(x.g); }
inline pari zk_to_ff(pari::tmp x, pari::tmp modpr)
{ return zk_to_ff(x.g, modpr.g); }
inline pari zkmodprinit(pari::tmp nf, pari::tmp pr)
{ return zkmodprinit(nf.g, pr.g); }
inline int RgV_isscalar(pari::tmp x)
{ return RgV_isscalar(x.g); }
inline pari algtobasis(pari::tmp nf, pari::tmp x)
{ return algtobasis(nf.g, x.g); }
inline pari arch_to_perm(pari::tmp arch)
{ return arch_to_perm(arch.g); }
inline pari basistoalg(pari::tmp nf, pari::tmp x)
{ return basistoalg(nf.g, x.g); }
inline pari dethnf(pari::tmp x)
{ return dethnf(x.g); }
inline pari dethnf_i(pari::tmp mat)
{ return dethnf_i(mat.g); }
inline pari element_div(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return element_div(nf.g, x.g, y.g); }
inline pari element_inv(pari::tmp nf, pari::tmp x)
{ return element_inv(nf.g, x.g); }
inline pari element_invmodideal(pari::tmp nf, pari::tmp x, pari::tmp ideal)
{ return element_invmodideal(nf.g, x.g, ideal.g); }
inline pari element_mul(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return element_mul(nf.g, x.g, y.g); }
inline pari element_muli(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return element_muli(nf.g, x.g, y.g); }
inline pari element_mulid(pari::tmp nf, pari::tmp x, long i)
{ return element_mulid(nf.g, x.g, i); }
inline pari element_pow(pari::tmp nf, pari::tmp x, pari::tmp k)
{ return element_pow(nf.g, x.g, k.g); }
inline pari element_powmodideal(pari::tmp nf, pari::tmp x, pari::tmp k, pari::tmp ideal)
{ return element_powmodideal(nf.g, x.g, k.g, ideal.g); }
inline pari element_powmodidele(pari::tmp nf, pari::tmp x, pari::tmp k, pari::tmp idele, pari::tmp structarch)
{ return element_powmodidele(nf.g, x.g, k.g, idele.g, structarch.g); }
inline pari element_sqr(pari::tmp nf, pari::tmp x)
{ return element_sqr(nf.g, x.g); }
inline pari element_sqri(pari::tmp nf, pari::tmp x)
{ return element_sqri(nf.g, x.g); }
inline long element_val(pari::tmp nf, pari::tmp x, pari::tmp vp)
{ return element_val(nf.g, x.g, vp.g); }
inline pari eltmul_get_table(pari::tmp nf, pari::tmp x)
{ return eltmul_get_table(nf.g, x.g); }
inline pari ideallist(pari::tmp nf, long bound)
{ return ideallist(nf.g, bound); }
inline pari ideallist0(pari::tmp nf, long bound, long flag)
{ return ideallist0(nf.g, bound, flag); }
inline pari ideallistarch(pari::tmp nf, pari::tmp list, pari::tmp arch)
{ return ideallistarch(nf.g, list.g, arch.g); }
inline pari ideallistunit(pari::tmp nf, long bound)
{ return ideallistunit(nf.g, bound); }
inline pari ideallistunitgen(pari::tmp nf, long bound)
{ return ideallistunitgen(nf.g, bound); }
inline pari ideallistzstar(pari::tmp nf, long bound)
{ return ideallistzstar(nf.g, bound); }
inline pari ideallistzstargen(pari::tmp nf, long bound)
{ return ideallistzstargen(nf.g, bound); }
inline pari idealstar0(pari::tmp nf, pari::tmp x, long flag)
{ return idealstar0(nf.g, x.g, flag); }
inline int isnfscalar(pari::tmp x)
{ return isnfscalar(x.g); }
inline pari lift_to_pol(pari::tmp x)
{ return lift_to_pol(x.g); }
inline pari lllreducemodmatrix(pari::tmp x, pari::tmp y)
{ return lllreducemodmatrix(x.g, y.g); }
inline pari matalgtobasis(pari::tmp nf, pari::tmp x)
{ return matalgtobasis(nf.g, x.g); }
inline pari matbasistoalg(pari::tmp nf, pari::tmp x)
{ return matbasistoalg(nf.g, x.g); }
inline pari nfdiveuc(pari::tmp nf, pari::tmp a, pari::tmp b)
{ return nfdiveuc(nf.g, a.g, b.g); }
inline pari nfdivrem(pari::tmp nf, pari::tmp a, pari::tmp b)
{ return nfdivrem(nf.g, a.g, b.g); }
inline pari nfmod(pari::tmp nf, pari::tmp a, pari::tmp b)
{ return nfmod(nf.g, a.g, b.g); }
inline pari nfreducemodideal(pari::tmp nf, pari::tmp x, pari::tmp ideal)
{ return nfreducemodideal(nf.g, x.g, ideal.g); }
inline pari nfreducemodidele(pari::tmp nf, pari::tmp g, pari::tmp idele, pari::tmp structarch)
{ return nfreducemodidele(nf.g, g.g, idele.g, structarch.g); }
inline pari reducemodinvertible(pari::tmp x, pari::tmp y)
{ return reducemodinvertible(x.g, y.g); }
inline pari reducemodmatrix(pari::tmp x, pari::tmp y)
{ return reducemodmatrix(x.g, y.g); }
inline pari rnfalgtobasis(pari::tmp rnf, pari::tmp x)
{ return rnfalgtobasis(rnf.g, x.g); }
inline pari rnfbasistoalg(pari::tmp rnf, pari::tmp x)
{ return rnfbasistoalg(rnf.g, x.g); }
inline pari set_sign_mod_idele(pari::tmp nf, pari::tmp x, pari::tmp y, pari::tmp idele, pari::tmp sarch)
{ return set_sign_mod_idele(nf.g, x.g, y.g, idele.g, sarch.g); }
inline pari vecmodii(pari::tmp a, pari::tmp b)
{ return vecmodii(a.g, b.g); }
inline pari zarchstar(pari::tmp nf, pari::tmp x, pari::tmp arch)
{ return zarchstar(nf.g, x.g, arch.g); }
inline pari zideallog(pari::tmp nf, pari::tmp x, pari::tmp bigideal)
{ return zideallog(nf.g, x.g, bigideal.g); }
inline pari zidealstar(pari::tmp nf, pari::tmp x)
{ return zidealstar(nf.g, x.g); }
inline pari zidealstarinit(pari::tmp nf, pari::tmp x)
{ return zidealstarinit(nf.g, x.g); }
inline pari Idealstar(pari::tmp nf, pari::tmp x, long flun)
{ return Idealstar(nf.g, x.g, flun); }
inline pari zidealstarinitgen(pari::tmp nf, pari::tmp x)
{ return zidealstarinitgen(nf.g, x.g); }
inline pari znlog(pari::tmp x, pari::tmp g)
{ return znlog(x.g, g.g); }
inline pari zsigne(pari::tmp nf, pari::tmp alpha, pari::tmp arch)
{ return zsigne(nf.g, alpha.g, arch.g); }
inline pari zsigns(pari::tmp nf, pari::tmp alpha)
{ return zsigns(nf.g, alpha.g); }
inline pari element_divmodpr(pari::tmp nf, pari::tmp x, pari::tmp y, pari::tmp modpr)
{ return element_divmodpr(nf.g, x.g, y.g, modpr.g); }
inline pari element_invmodpr(pari::tmp nf, pari::tmp y, pari::tmp modpr)
{ return element_invmodpr(nf.g, y.g, modpr.g); }
inline pari element_mulmodpr(pari::tmp nf, pari::tmp x, pari::tmp y, pari::tmp modpr)
{ return element_mulmodpr(nf.g, x.g, y.g, modpr.g); }
inline pari element_mulvec(pari::tmp nf, pari::tmp x, pari::tmp v)
{ return element_mulvec(nf.g, x.g, v.g); }
inline pari element_powmodpr(pari::tmp nf, pari::tmp x, pari::tmp k, pari::tmp modpr)
{ return element_powmodpr(nf.g, x.g, k.g, modpr.g); }
inline pari element_reduce(pari::tmp nf, pari::tmp x, pari::tmp ideal)
{ return element_reduce(nf.g, x.g, ideal.g); }
inline pari ideal_two_elt(pari::tmp nf, pari::tmp ix)
{ return ideal_two_elt(nf.g, ix.g); }
inline pari ideal_two_elt0(pari::tmp nf, pari::tmp ix, pari::tmp a)
{ return ideal_two_elt0(nf.g, ix.g, a.g); }
inline pari ideal_two_elt2(pari::tmp nf, pari::tmp x, pari::tmp a)
{ return ideal_two_elt2(nf.g, x.g, a.g); }
inline pari idealadd(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return idealadd(nf.g, x.g, y.g); }
inline pari idealaddmultoone(pari::tmp nf, pari::tmp list)
{ return idealaddmultoone(nf.g, list.g); }
inline pari idealaddtoone(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return idealaddtoone(nf.g, x.g, y.g); }
inline pari idealaddtoone0(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return idealaddtoone0(nf.g, x.g, y.g); }
inline pari idealappr(pari::tmp nf, pari::tmp x)
{ return idealappr(nf.g, x.g); }
inline pari idealappr0(pari::tmp nf, pari::tmp x, long fl)
{ return idealappr0(nf.g, x.g, fl); }
inline pari idealapprfact(pari::tmp nf, pari::tmp x)
{ return idealapprfact(nf.g, x.g); }
inline pari idealchinese(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return idealchinese(nf.g, x.g, y.g); }
inline pari idealcoprime(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return idealcoprime(nf.g, x.g, y.g); }
inline pari idealdiv(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return idealdiv(nf.g, x.g, y.g); }
inline pari idealdiv0(pari::tmp nf, pari::tmp x, pari::tmp y, long flag)
{ return idealdiv0(nf.g, x.g, y.g, flag); }
inline pari idealdivexact(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return idealdivexact(nf.g, x.g, y.g); }
inline pari idealdivpowprime(pari::tmp nf, pari::tmp x, pari::tmp vp, pari::tmp n)
{ return idealdivpowprime(nf.g, x.g, vp.g, n.g); }
inline pari idealmulpowprime(pari::tmp nf, pari::tmp x, pari::tmp vp, pari::tmp n)
{ return idealmulpowprime(nf.g, x.g, vp.g, n.g); }
inline pari idealfactor(pari::tmp nf, pari::tmp x)
{ return idealfactor(nf.g, x.g); }
inline pari idealhermite(pari::tmp nf, pari::tmp x)
{ return idealhermite(nf.g, x.g); }
inline pari idealhnf0(pari::tmp nf, pari::tmp a, pari::tmp b)
{ return idealhnf0(nf.g, a.g, b.g); }
inline pari idealintersect(pari::tmp nf, pari::tmp x, pari::tmp y)
{ return idealintersect(nf.g, x.g, y.g); }
inline pari idealinv(pari::tmp nf, pari::tmp ix)
{ return idealinv(nf.g, ix.g); }
inline pari ideallllred(pari::tmp nf, pari::tmp ix, pari::tmp vdir, long prc = prec)
{ return ideallllred(nf.g, ix.g, vdir.g, prc); }
inline pari ideallllred_elt(pari::tmp nf, pari::tmp I, pari::tmp vdir)
{ return ideallllred_elt(nf.g, I.g, vdir.g); }
inline pari idealmul(pari::tmp nf, pari::tmp ix, pari::tmp iy)
{ return idealmul(nf.g, ix.g, iy.g); }
inline pari idealmul0(pari::tmp nf, pari::tmp ix, pari::tmp iy, long flag, long prc = prec)
{ return idealmul0(nf.g, ix.g, iy.g, flag, prc); }
inline pari idealmulh(pari::tmp nf, pari::tmp ix, pari::tmp iy)
{ return idealmulh(nf.g, ix.g, iy.g); }
inline pari idealmulprime(pari::tmp nf, pari::tmp ix, pari::tmp vp)
{ return idealmulprime(nf.g, ix.g, vp.g); }
inline pari idealmulred(pari::tmp nf, pari::tmp ix, pari::tmp iy, long prc = prec)
{ return idealmulred(nf.g, ix.g, iy.g, prc); }
inline pari idealnorm(pari::tmp nf, pari::tmp x)
{ return idealnorm(nf.g, x.g); }
inline pari idealpow(pari::tmp nf, pari::tmp ix, pari::tmp n)
{ return idealpow(nf.g, ix.g, n.g); }
inline pari idealpow0(pari::tmp nf, pari::tmp ix, pari::tmp n, long flag, long prc = prec)
{ return idealpow0(nf.g, ix.g, n.g, flag, prc); }
inline pari idealpowred(pari::tmp nf, pari::tmp ix, pari::tmp n, long prc = prec)
{ return idealpowred(nf.g, ix.g, n.g, prc); }
inline pari idealpows(pari::tmp nf, pari::tmp ideal, long iexp)
{ return idealpows(nf.g, ideal.g, iexp); }
inline pari idealprodprime(pari::tmp nf, pari::tmp L)
{ return idealprodprime(nf.g, L.g); }
inline pari idealred_elt(pari::tmp nf, pari::tmp I)
{ return idealred_elt(nf.g, I.g); }
inline long idealval(pari::tmp nf, pari::tmp ix, pari::tmp vp)
{ return idealval(nf.g, ix.g, vp.g); }
inline pari ideleaddone(pari::tmp nf, pari::tmp x, pari::tmp idele)
{ return ideleaddone(nf.g, x.g, idele.g); }
inline int isidentity(pari::tmp x)
{ return isidentity(x.g); }
inline long isideal(pari::tmp nf, pari::tmp x)
{ return isideal(nf.g, x.g); }
inline pari minideal(pari::tmp nf, pari::tmp ix, pari::tmp vdir, long prc = prec)
{ return minideal(nf.g, ix.g, vdir.g, prc); }
inline pari mul_content(pari::tmp cx, pari::tmp cy)
{ return mul_content(cx.g, cy.g); }
inline pari nfdetint(pari::tmp nf, pari::tmp pseudo)
{ return nfdetint(nf.g, pseudo.g); }
inline pari nfhermite(pari::tmp nf, pari::tmp x)
{ return nfhermite(nf.g, x.g); }
inline pari nfhermitemod(pari::tmp nf, pari::tmp x, pari::tmp detmat)
{ return nfhermitemod(nf.g, x.g, detmat.g); }
inline pari nfkermodpr(pari::tmp nf, pari::tmp x, pari::tmp modpr)
{ return nfkermodpr(nf.g, x.g, modpr.g); }
inline pari nfsmith(pari::tmp nf, pari::tmp x)
{ return nfsmith(nf.g, x.g); }
inline pari nfsolvemodpr(pari::tmp nf, pari::tmp a, pari::tmp b, pari::tmp modpr)
{ return nfsolvemodpr(nf.g, a.g, b.g, modpr.g); }
inline pari prime_to_ideal(pari::tmp nf, pari::tmp vp)
{ return prime_to_ideal(nf.g, vp.g); }
inline pari principalideal(pari::tmp nf, pari::tmp a)
{ return principalideal(nf.g, a.g); }
inline pari principalidele(pari::tmp nf, pari::tmp a, long prc = prec)
{ return principalidele(nf.g, a.g, prc); }
inline pari vecdiv(pari::tmp x, pari::tmp y)
{ return vecdiv(x.g, y.g); }
inline pari vecinv(pari::tmp x)
{ return vecinv(x.g); }
inline pari vecmul(pari::tmp x, pari::tmp y)
{ return vecmul(x.g, y.g); }
inline pari vecpow(pari::tmp x, pari::tmp n)
{ return vecpow(x.g, n.g); }
inline pari rnfelementabstorel(pari::tmp rnf, pari::tmp x)
{ return rnfelementabstorel(rnf.g, x.g); }
inline pari rnfelementdown(pari::tmp rnf, pari::tmp x)
{ return rnfelementdown(rnf.g, x.g); }
inline pari rnfelementreltoabs(pari::tmp rnf, pari::tmp x)
{ return rnfelementreltoabs(rnf.g, x.g); }
inline pari rnfelementup(pari::tmp rnf, pari::tmp x)
{ return rnfelementup(rnf.g, x.g); }
inline pari rnfidealabstorel(pari::tmp rnf, pari::tmp x)
{ return rnfidealabstorel(rnf.g, x.g); }
inline pari rnfidealdown(pari::tmp rnf, pari::tmp x)
{ return rnfidealdown(rnf.g, x.g); }
inline pari rnfidealhermite(pari::tmp rnf, pari::tmp x)
{ return rnfidealhermite(rnf.g, x.g); }
inline pari rnfidealmul(pari::tmp rnf, pari::tmp x, pari::tmp y)
{ return rnfidealmul(rnf.g, x.g, y.g); }
inline pari rnfidealnormabs(pari::tmp rnf, pari::tmp x)
{ return rnfidealnormabs(rnf.g, x.g); }
inline pari rnfidealnormrel(pari::tmp rnf, pari::tmp x)
{ return rnfidealnormrel(rnf.g, x.g); }
inline pari rnfidealreltoabs(pari::tmp rnf, pari::tmp x)
{ return rnfidealreltoabs(rnf.g, x.g); }
inline pari rnfidealtwoelement(pari::tmp rnf, pari::tmp x)
{ return rnfidealtwoelement(rnf.g, x.g); }
inline pari rnfidealup(pari::tmp rnf, pari::tmp x)
{ return rnfidealup(rnf.g, x.g); }
inline pari rnfinitalg(pari::tmp nf, pari::tmp pol, long prc = prec)
{ return rnfinitalg(nf.g, pol.g, prc); }
inline pari ZM_zc_mul(pari::tmp x, pari::tmp y)
{ return ZM_zc_mul(x.g, y.g); }
inline pari ZM_zm_mul(pari::tmp x, pari::tmp y)
{ return ZM_zm_mul(x.g, y.g); }
inline pari algdep(pari::tmp x, long n, long prc = prec)
{ return algdep(x.g, n, prc); }
inline pari algdep0(pari::tmp x, long n, long bit, long prc = prec)
{ return algdep0(x.g, n, bit, prc); }
inline pari algdep2(pari::tmp x, long n, long bit)
{ return algdep2(x.g, n, bit); }
inline pari gram_matrix(pari::tmp M)
{ return gram_matrix(M.g); }
inline pari kerint(pari::tmp x)
{ return kerint(x.g); }
inline pari kerint1(pari::tmp x)
{ return kerint1(x.g); }
inline pari lindep(pari::tmp x, long prc = prec)
{ return lindep(x.g, prc); }
inline pari lindep0(pari::tmp x, long flag, long prc = prec)
{ return lindep0(x.g, flag, prc); }
inline pari lindep2(pari::tmp x, long bit)
{ return lindep2(x.g, bit); }
inline pari lll(pari::tmp x, long prc = prec)
{ return lll(x.g, prc); }
inline pari lllgen(pari::tmp x)
{ return lllgen(x.g); }
inline pari lllgram(pari::tmp x, long prc = prec)
{ return lllgram(x.g, prc); }
inline pari lllgramgen(pari::tmp x)
{ return lllgramgen(x.g); }
inline pari lllgramint(pari::tmp x)
{ return lllgramint(x.g); }
inline pari lllgramkerim(pari::tmp x)
{ return lllgramkerim(x.g); }
inline pari lllgramkerimgen(pari::tmp x)
{ return lllgramkerimgen(x.g); }
inline pari lllint(pari::tmp x)
{ return lllint(x.g); }
inline pari lllint_ip(pari::tmp x, long alpha)
{ return lllint_ip(x.g, alpha); }
inline pari lllintpartial(pari::tmp mat)
{ return lllintpartial(mat.g); }
inline pari lllintpartial_ip(pari::tmp mat)
{ return lllintpartial_ip(mat.g); }
inline pari lllkerim(pari::tmp x)
{ return lllkerim(x.g); }
inline pari lllkerimgen(pari::tmp x)
{ return lllkerimgen(x.g); }
inline pari matkerint0(pari::tmp x, long flag)
{ return matkerint0(x.g, flag); }
inline pari minim(pari::tmp a, pari::tmp borne, pari::tmp stockmax)
{ return minim(a.g, borne.g, stockmax.g); }
inline pari qfrep0(pari::tmp a, pari::tmp borne, long flag)
{ return qfrep0(a.g, borne.g, flag); }
inline pari qfminim0(pari::tmp a, pari::tmp borne, pari::tmp stockmax, long flag, long prc = prec)
{ return qfminim0(a.g, borne.g, stockmax.g, flag, prc); }
inline pari minim2(pari::tmp a, pari::tmp borne, pari::tmp stockmax)
{ return minim2(a.g, borne.g, stockmax.g); }
inline pari perf(pari::tmp a)
{ return perf(a.g); }
inline pari qflll0(pari::tmp x, long flag, long prc = prec)
{ return qflll0(x.g, flag, prc); }
inline pari qflllgram0(pari::tmp x, long flag, long prc = prec)
{ return qflllgram0(x.g, flag, prc); }
inline pari binomial(pari::tmp x, long k)
{ return binomial(x.g, k); }
inline int cmp_prime_ideal(pari::tmp x, pari::tmp y)
{ return cmp_prime_ideal(x.g, y.g); }
inline int cmp_prime_over_p(pari::tmp x, pari::tmp y)
{ return cmp_prime_over_p(x.g, y.g); }
inline int cmp_vecint(pari::tmp x, pari::tmp y)
{ return cmp_vecint(x.g, y.g); }
inline pari convol(pari::tmp x, pari::tmp y)
{ return convol(x.g, y.g); }
inline pari p_cyclo(long n, long v)
{ return cyclo(n, v); }
inline pari dirdiv(pari::tmp x, pari::tmp y)
{ return dirdiv(x.g, y.g); }
inline pari dirmul(pari::tmp x, pari::tmp y)
{ return dirmul(x.g, y.g); }
inline pari genrand(pari::tmp N)
{ return genrand(N.g); }
inline pari gprec(pari::tmp x, long l)
{ return gprec(x.g, l); }
inline pari gprec_wtrunc(pari::tmp x, long pr)
{ return gprec_wtrunc(x.g, pr); }
inline pari gprec_w(pari::tmp x, long pr)
{ return gprec_w(x.g, pr); }
inline pari gtoset(pari::tmp x)
{ return gtoset(x.g); }
inline pari indexlexsort(pari::tmp x)
{ return indexlexsort(x.g); }
inline pari indexsort(pari::tmp x)
{ return indexsort(x.g); }
inline pari laplace(pari::tmp x)
{ return laplace(x.g); }
inline pari p_legendre(long n, long v)
{ return legendre(n, v); }
inline pari lexsort(pari::tmp x)
{ return lexsort(x.g); }
inline pari p_mathilbert(long n)
{ return mathilbert(n); }
inline pari matqpascal(long n, pari::tmp q)
{ return matqpascal(n, q.g); }
inline pari modreverse_i(pari::tmp a, pari::tmp T)
{ return modreverse_i(a.g, T.g); }
inline pari numtoperm(long n, pari::tmp x)
{ return numtoperm(n, x.g); }
inline int pari_compare_lg(pari::tmp a, pari::tmp b)
{ return pari_compare_lg(a.g, b.g); }
inline pari permtonum(pari::tmp x)
{ return permtonum(x.g); }
inline pari polrecip(pari::tmp x)
{ return polrecip(x.g); }
inline pari polymodrecip(pari::tmp x)
{ return polymodrecip(x.g); }
inline pari roots_to_pol(pari::tmp a, long v)
{ return roots_to_pol(a.g, v); }
inline pari setintersect(pari::tmp x, pari::tmp y)
{ return setintersect(x.g, y.g); }
inline long setisset(pari::tmp x)
{ return setisset(x.g); }
inline pari setminus(pari::tmp x, pari::tmp y)
{ return setminus(x.g, y.g); }
inline long setsearch(pari::tmp x, pari::tmp y, long flag)
{ return setsearch(x.g, y.g, flag); }
inline pari setunion(pari::tmp x, pari::tmp y)
{ return setunion(x.g, y.g); }
inline pari sindexlexsort(pari::tmp x)
{ return sindexlexsort(x.g); }
inline pari sindexsort(pari::tmp x)
{ return sindexsort(x.g); }
inline pari sort(pari::tmp x)
{ return sort(x.g); }
inline pari p_tchebi(long n, long v)
{ return tchebi(n, v); }
inline pari p_vecbinome(long n)
{ return vecbinome(n); }
inline pari vecsort(pari::tmp x, pari::tmp k)
{ return vecsort(x.g, k.g); }
inline pari vecsort0(pari::tmp x, pari::tmp k, long flag)
{ return vecsort0(x.g, k.g, flag); }
inline pari buchimag(pari::tmp D, pari::tmp gcbach, pari::tmp gcbach2, pari::tmp gCO)
{ return buchimag(D.g, gcbach.g, gcbach2.g, gCO.g); }
inline pari buchreal(pari::tmp D, pari::tmp gsens, pari::tmp gcbach, pari::tmp gcbach2, pari::tmp gRELSUP, long prc = prec)
{ return buchreal(D.g, gsens.g, gcbach.g, gcbach2.g, gRELSUP.g, prc); }
inline pari quadclassunit0(pari::tmp x, long flag, pari::tmp data, long prc = prec)
{ return quadclassunit0(x.g, flag, data.g, prc); }
inline pari quadhilbert(pari::tmp D, pari::tmp flag, long prc = prec)
{ return quadhilbert(D.g, flag.g, prc); }
inline pari quadray(pari::tmp bnf, pari::tmp f, pari::tmp flag, long prc = prec)
{ return quadray(bnf.g, f.g, flag.g, prc); }
inline pari bnfclassunit0(pari::tmp P, long flag, pari::tmp data, long prc = prec)
{ return bnfclassunit0(P.g, flag, data.g, prc); }
inline pari bnfinit0(pari::tmp P, long flag, pari::tmp data, long prc = prec)
{ return bnfinit0(P.g, flag, data.g, prc); }
inline pari bnfmake(pari::tmp sbnf, long prc = prec)
{ return bnfmake(sbnf.g, prc); }
inline pari bnfnewprec(pari::tmp nf, long prc = prec)
{ return bnfnewprec(nf.g, prc); }
inline pari bnrnewprec(pari::tmp bnr, long prc = prec)
{ return bnrnewprec(bnr.g, prc); }
inline pari buchall(pari::tmp P, double bach, double bach2, long nbrelpid, long flun, long prc = prec)
{ return buchall(P.g, bach, bach2, nbrelpid, flun, prc); }
inline pari buchfu(pari::tmp bignf)
{ return buchfu(bignf.g); }
inline pari classgrouponly(pari::tmp P, pari::tmp data, long prc = prec)
{ return classgrouponly(P.g, data.g, prc); }
inline pari isprincipal(pari::tmp bignf, pari::tmp x)
{ return isprincipal(bignf.g, x.g); }
inline pari isprincipalall(pari::tmp bignf, pari::tmp x, long flall)
{ return isprincipalall(bignf.g, x.g, flall); }
inline pari isprincipalfact(pari::tmp bnf, pari::tmp P, pari::tmp e, pari::tmp C, long flag)
{ return isprincipalfact(bnf.g, P.g, e.g, C.g, flag); }
inline pari isprincipalforce(pari::tmp bignf, pari::tmp x)
{ return isprincipalforce(bignf.g, x.g); }
inline pari isprincipalgen(pari::tmp bignf, pari::tmp x)
{ return isprincipalgen(bignf.g, x.g); }
inline pari isprincipalgenforce(pari::tmp bignf, pari::tmp x)
{ return isprincipalgenforce(bignf.g, x.g); }
inline pari isunit(pari::tmp bignf, pari::tmp x)
{ return isunit(bignf.g, x.g); }
inline pari regulator(pari::tmp P, pari::tmp data, long prc = prec)
{ return regulator(P.g, data.g, prc); }
inline pari signunits(pari::tmp bignf)
{ return signunits(bignf.g); }
inline pari smallbuchinit(pari::tmp pol, double bach, double bach2, long nbrelpid, long prc = prec)
{ return smallbuchinit(pol.g, bach, bach2, nbrelpid, prc); }
inline pari zsignunits(pari::tmp bnf, pari::tmp archp, int add_zu)
{ return zsignunits(bnf.g, archp.g, add_zu); }
inline pari bnrclass0(pari::tmp bignf, pari::tmp ideal, long flag)
{ return bnrclass0(bignf.g, ideal.g, flag); }
inline pari bnrclassno(pari::tmp bignf, pari::tmp ideal)
{ return bnrclassno(bignf.g, ideal.g); }
inline pari bnrclassnolist(pari::tmp bnf, pari::tmp listes)
{ return bnrclassnolist(bnf.g, listes.g); }
inline pari bnrconductor(pari::tmp arg0, pari::tmp arg1, pari::tmp arg2, pari::tmp flag)
{ return bnrconductor(arg0.g, arg1.g, arg2.g, flag.g); }
inline pari bnrconductorofchar(pari::tmp bnr, pari::tmp chi)
{ return bnrconductorofchar(bnr.g, chi.g); }
inline pari bnrdisc0(pari::tmp arg0, pari::tmp arg1, pari::tmp arg2, long flag)
{ return bnrdisc0(arg0.g, arg1.g, arg2.g, flag); }
inline pari bnrdisclist0(pari::tmp bnf, pari::tmp borne, pari::tmp arch)
{ return bnrdisclist0(bnf.g, borne.g, arch.g); }
inline pari bnrinit0(pari::tmp bignf, pari::tmp ideal, long flag)
{ return bnrinit0(bignf.g, ideal.g, flag); }
inline long bnrisconductor(pari::tmp arg0, pari::tmp arg1, pari::tmp arg2)
{ return bnrisconductor(arg0.g, arg1.g, arg2.g); }
inline pari bnrisprincipal(pari::tmp bnf, pari::tmp x, long flag)
{ return bnrisprincipal(bnf.g, x.g, flag); }
inline pari buchnarrow(pari::tmp bignf)
{ return buchnarrow(bignf.g); }
inline pari buchray(pari::tmp bignf, pari::tmp ideal)
{ return buchray(bignf.g, ideal.g); }
inline pari buchrayinit(pari::tmp bignf, pari::tmp ideal)
{ return buchrayinit(bignf.g, ideal.g); }
inline pari buchrayinitgen(pari::tmp bignf, pari::tmp ideal)
{ return buchrayinitgen(bignf.g, ideal.g); }
inline long certifybuchall(pari::tmp bnf)
{ return certifybuchall(bnf.g); }
inline pari conductor(pari::tmp bnr, pari::tmp subgroup, long all)
{ return conductor(bnr.g, subgroup.g, all); }
inline pari decodemodule(pari::tmp nf, pari::tmp fa)
{ return decodemodule(nf.g, fa.g); }
inline pari discrayabs(pari::tmp bnr, pari::tmp subgroup)
{ return discrayabs(bnr.g, subgroup.g); }
inline pari discrayabscond(pari::tmp bnr, pari::tmp subgroup)
{ return discrayabscond(bnr.g, subgroup.g); }
inline pari discrayabslist(pari::tmp bnf, pari::tmp listes)
{ return discrayabslist(bnf.g, listes.g); }
inline pari discrayabslistarch(pari::tmp bnf, pari::tmp arch, long bound)
{ return discrayabslistarch(bnf.g, arch.g, bound); }
inline pari discrayabslistlong(pari::tmp bnf, long bound)
{ return discrayabslistlong(bnf.g, bound); }
inline pari discrayrel(pari::tmp bnr, pari::tmp subgroup)
{ return discrayrel(bnr.g, subgroup.g); }
inline pari discrayrelcond(pari::tmp bnr, pari::tmp subgroup)
{ return discrayrelcond(bnr.g, subgroup.g); }
inline pari idealmodidele(pari::tmp bnr, pari::tmp x)
{ return idealmodidele(bnr.g, x.g); }
inline long isinvector(pari::tmp v, pari::tmp x)
{ return isinvector(v.g, x.g); }
inline pari isprincipalray(pari::tmp bnf, pari::tmp x)
{ return isprincipalray(bnf.g, x.g); }
inline pari isprincipalraygen(pari::tmp bnf, pari::tmp x)
{ return isprincipalraygen(bnf.g, x.g); }
inline pari quick_isprincipalgen(pari::tmp bnf, pari::tmp x)
{ return quick_isprincipalgen(bnf.g, x.g); }
inline pari rnfconductor(pari::tmp bnf, pari::tmp polrel, long flag)
{ return rnfconductor(bnf.g, polrel.g, flag); }
inline pari rnfnormgroup(pari::tmp bnr, pari::tmp polrel)
{ return rnfnormgroup(bnr.g, polrel.g); }
inline pari subgrouplist0(pari::tmp bnr, pari::tmp indexbound, long all)
{ return subgrouplist0(bnr.g, indexbound.g, all); }
inline pari bnfisnorm(pari::tmp bnf, pari::tmp x, long flag, long PREC)
{ return bnfisnorm(bnf.g, x.g, flag, PREC); }
inline pari rnfisnorm(pari::tmp S, pari::tmp x, long flag)
{ return rnfisnorm(S.g, x.g, flag); }
inline pari rnfisnorminit(pari::tmp bnf, pari::tmp relpol, int galois)
{ return rnfisnorminit(bnf.g, relpol.g, galois); }
inline pari bnfissunit(pari::tmp bnf, pari::tmp suni, pari::tmp x)
{ return bnfissunit(bnf.g, suni.g, x.g); }
inline pari bnfsunit(pari::tmp bnf, pari::tmp s, long PREC)
{ return bnfsunit(bnf.g, s.g, PREC); }
inline long nfhilbert(pari::tmp bnf, pari::tmp a, pari::tmp b)
{ return nfhilbert(bnf.g, a.g, b.g); }
inline long nfhilbert0(pari::tmp bnf, pari::tmp a, pari::tmp b, pari::tmp p)
{ return nfhilbert0(bnf.g, a.g, b.g, p.g); }
inline long nfhilbertp(pari::tmp bnf, pari::tmp a, pari::tmp b, pari::tmp p)
{ return nfhilbertp(bnf.g, a.g, b.g, p.g); }
inline long qpsoluble(pari::tmp pol, pari::tmp p)
{ return qpsoluble(pol.g, p.g); }
inline long qpsolublenf(pari::tmp bnf, pari::tmp pol, pari::tmp p)
{ return qpsolublenf(bnf.g, pol.g, p.g); }
inline long zpsoluble(pari::tmp pol, pari::tmp p)
{ return zpsoluble(pol.g, p.g); }
inline long zpsolublenf(pari::tmp bnf, pari::tmp pol, pari::tmp p)
{ return zpsolublenf(bnf.g, pol.g, p.g); }
inline pari p_default0(char * a, char * b, long flag)
{ return default0(a, b, flag); }
inline pari p_gp_default(char * a, char * b)
{ return gp_default(a, b); }
inline pari p_ellcondfile(long f)
{ return ellcondfile(f); }
inline pari p_ellcondlist(long f)
{ return ellcondlist(f); }
inline pari ellconvertname(pari::tmp s)
{ return ellconvertname(s.g); }
inline pari ellgenerators(pari::tmp E)
{ return ellgenerators(E.g); }
inline pari ellidentify(pari::tmp E)
{ return ellidentify(E.g); }
inline pari ellsearch(pari::tmp A)
{ return ellsearch(A.g); }
inline pari ellsearchcurve(pari::tmp name)
{ return ellsearchcurve(name.g); }
inline pari addell(pari::tmp e, pari::tmp z1, pari::tmp z2)
{ return addell(e.g, z1.g, z2.g); }
inline pari akell(pari::tmp e, pari::tmp n)
{ return akell(e.g, n.g); }
inline pari anell(pari::tmp e, long n)
{ return anell(e.g, n); }
inline pari apell(pari::tmp e, pari::tmp p)
{ return apell(e.g, p.g); }
inline pari apell2(pari::tmp e, pari::tmp p)
{ return apell2(e.g, p.g); }
inline pari bilhell(pari::tmp e, pari::tmp z1, pari::tmp z2, long prc = prec)
{ return bilhell(e.g, z1.g, z2.g, prc); }
inline void checkbell(pari::tmp e)
{ checkbell(e.g); }
inline void checkell(pari::tmp e)
{ checkell(e.g); }
inline void checksell(pari::tmp e)
{ checksell(e.g); }
inline pari coordch(pari::tmp e, pari::tmp ch)
{ return coordch(e.g, ch.g); }
inline pari ellap0(pari::tmp e, pari::tmp p, long flag)
{ return ellap0(e.g, p.g, flag); }
inline pari elleisnum(pari::tmp om, long k, long flag, long prc = prec)
{ return elleisnum(om.g, k, flag, prc); }
inline pari elleta(pari::tmp om, long prc = prec)
{ return elleta(om.g, prc); }
inline pari ellglobalred(pari::tmp e1)
{ return ellglobalred(e1.g); }
inline pari ellheight0(pari::tmp e, pari::tmp a, long flag, long prc = prec)
{ return ellheight0(e.g, a.g, flag, prc); }
inline pari ellinit0(pari::tmp x, long flag, long prc = prec)
{ return ellinit0(x.g, flag, prc); }
inline pari ellisoncurve(pari::tmp e, pari::tmp z)
{ return ellisoncurve(e.g, z.g); }
inline pari elllseries(pari::tmp e, pari::tmp s, pari::tmp A, long prc = prec)
{ return elllseries(e.g, s.g, A.g, prc); }
inline pari elllocalred(pari::tmp e, pari::tmp p1)
{ return elllocalred(e.g, p1.g); }
inline long ellrootno(pari::tmp e, pari::tmp p)
{ return ellrootno(e.g, p.g); }
inline pari ellsigma(pari::tmp om, pari::tmp z, long flag, long prc = prec)
{ return ellsigma(om.g, z.g, flag, prc); }
inline pari elltaniyama(pari::tmp e, long prc = prec)
{ return elltaniyama(e.g, prc); }
inline pari elltors0(pari::tmp e, long flag)
{ return elltors0(e.g, flag); }
inline pari ellwp0(pari::tmp e, pari::tmp z, long flag, long prec, long PREC)
{ return ellwp0(e.g, z.g, flag, prec, PREC); }
inline pari ellzeta(pari::tmp om, pari::tmp z, long prc = prec)
{ return ellzeta(om.g, z.g, prc); }
inline pari ghell(pari::tmp e, pari::tmp a, long prc = prec)
{ return ghell(e.g, a.g, prc); }
inline pari ghell2(pari::tmp e, pari::tmp a, long prc = prec)
{ return ghell2(e.g, a.g, prc); }
inline pari initell(pari::tmp x, long prc = prec)
{ return initell(x.g, prc); }
inline pari mathell(pari::tmp e, pari::tmp x, long prc = prec)
{ return mathell(e.g, x.g, prc); }
inline int oncurve(pari::tmp e, pari::tmp z)
{ return oncurve(e.g, z.g); }
inline pari ordell(pari::tmp e, pari::tmp x, long prc = prec)
{ return ordell(e.g, x.g, prc); }
inline pari orderell(pari::tmp e, pari::tmp p)
{ return orderell(e.g, p.g); }
inline pari pointch(pari::tmp x, pari::tmp ch)
{ return pointch(x.g, ch.g); }
inline pari pointchinv(pari::tmp x, pari::tmp ch)
{ return pointchinv(x.g, ch.g); }
inline pari pointell(pari::tmp e, pari::tmp z, long prc = prec)
{ return pointell(e.g, z.g, prc); }
inline pari powell(pari::tmp e, pari::tmp z, pari::tmp n)
{ return powell(e.g, z.g, n.g); }
inline pari smallinitell(pari::tmp x)
{ return smallinitell(x.g); }
inline pari subell(pari::tmp e, pari::tmp z1, pari::tmp z2)
{ return subell(e.g, z1.g, z2.g); }
inline pari torsell(pari::tmp e)
{ return torsell(e.g); }
inline pari weipell(pari::tmp e, long precdl)
{ return weipell(e.g, precdl); }
inline pari zell(pari::tmp e, pari::tmp z, long prc = prec)
{ return zell(e.g, z.g, prc); }
inline pari GENtocanonicalstr(pari::tmp x)
{ return GENtocanonicalstr(x.g); }
inline pari GENtoGENstr(pari::tmp x)
{ return GENtoGENstr(x.g); }
inline pari Str(pari::tmp g)
{ return Str(g.g); }
inline pari Strchr(pari::tmp g)
{ return Strchr(g.g); }
inline pari Strexpand(pari::tmp g)
{ return Strexpand(g.g); }
inline pari Strtex(pari::tmp g)
{ return Strtex(g.g); }
inline void brute(pari::tmp g, char format, long dec)
{ brute(g.g, format, dec); }
inline void bruteall(pari::tmp g, char f, long d, long sp)
{ bruteall(g.g, f, d, sp); }
inline void bruterr(pari::tmp x, char format, long dec)
{ bruterr(x.g, format, dec); }
inline void error0(pari::tmp g)
{ error0(g.g); }
inline pari p_gp_read_file(char * s)
{ return gp_read_file(s); }
inline pari p_gp_read_stream(FILE * f)
{ return gp_read_stream(f); }
inline pari p_gp_readvec_file(char * s)
{ return gp_readvec_file(s); }
inline pari p_gp_readvec_stream(FILE * f)
{ return gp_readvec_stream(f); }
inline void matbrute(pari::tmp g, char format, long dec)
{ matbrute(g.g, format, dec); }
inline void outbeaut(pari::tmp x)
{ outbeaut(x.g); }
inline void outbeauterr(pari::tmp x)
{ outbeauterr(x.g); }
inline void outbrute(pari::tmp x)
{ outbrute(x.g); }
inline void outerr(pari::tmp x)
{ outerr(x.g); }
inline void outmat(pari::tmp x)
{ outmat(x.g); }
inline void output(pari::tmp x)
{ output(x.g); }
inline void outsor(pari::tmp x)
{ outsor(x.g); }
inline void outtex(pari::tmp x)
{ outtex(x.g); }
inline void print1(pari::tmp g)
{ print1(g.g); }
inline void printp(pari::tmp g)
{ printp(g.g); }
inline void printp1(pari::tmp g)
{ printp1(g.g); }
inline void printtex(pari::tmp g)
{ printtex(g.g); }
inline void sor(pari::tmp g, char fo, long dd, long chmp)
{ sor(g.g, fo, dd, chmp); }
inline void texe(pari::tmp g, char format, long dec)
{ texe(g.g, format, dec); }
inline void voir(pari::tmp x, long nb)
{ voir(x.g, nb); }
inline void write0(const char * s, pari::tmp g)
{ write0(s, g.g); }
inline void write1(const char * s, pari::tmp g)
{ write1(s, g.g); }
inline void writebin(char * name, pari::tmp x)
{ writebin(name, x.g); }
inline void writetex(const char * s, pari::tmp g)
{ writetex(s, g.g); }
inline pari checkgal(pari::tmp gal)
{ return checkgal(gal.g); }
inline pari galois_group(pari::tmp gal)
{ return galois_group(gal.g); }
inline pari galoisconj(pari::tmp nf)
{ return galoisconj(nf.g); }
inline pari galoisconj0(pari::tmp nf, long flag, pari::tmp d, long prc = prec)
{ return galoisconj0(nf.g, flag, d.g, prc); }
inline pari galoisconj2(pari::tmp x, long nbmax, long prc = prec)
{ return galoisconj2(x.g, nbmax, prc); }
inline pari galoisconj4(pari::tmp T, pari::tmp den, long flag)
{ return galoisconj4(T.g, den.g, flag); }
inline pari galoisexport(pari::tmp gal, long format)
{ return galoisexport(gal.g, format); }
inline pari galoisfixedfield(pari::tmp gal, pari::tmp v, long flag, long y)
{ return galoisfixedfield(gal.g, v.g, flag, y); }
inline pari galoisidentify(pari::tmp gal)
{ return galoisidentify(gal.g); }
inline pari galoisinit(pari::tmp nf, pari::tmp den)
{ return galoisinit(nf.g, den.g); }
inline pari galoisisabelian(pari::tmp gal, long flag)
{ return galoisisabelian(gal.g, flag); }
inline pari galoispermtopol(pari::tmp gal, pari::tmp perm)
{ return galoispermtopol(gal.g, perm.g); }
inline pari galoissubgroups(pari::tmp G)
{ return galoissubgroups(G.g); }
inline pari galoissubfields(pari::tmp G, long flag, long v)
{ return galoissubfields(G.g, flag, v); }
inline long numberofconjugates(pari::tmp T, long pdepart)
{ return numberofconjugates(T.g, pdepart); }
inline pari vandermondeinverse(pari::tmp L, pari::tmp T, pari::tmp den, pari::tmp prep)
{ return vandermondeinverse(L.g, T.g, den.g, prep.g); }
inline pari gadd(pari::tmp x, pari::tmp y)
{ return gadd(x.g, y.g); }
inline pari gaddsg(long x, pari::tmp y)
{ return gaddsg(x, y.g); }
inline pari gdiv(pari::tmp x, pari::tmp y)
{ return gdiv(x.g, y.g); }
inline pari gdivgs(pari::tmp x, long s)
{ return gdivgs(x.g, s); }
inline pari gmul(pari::tmp x, pari::tmp y)
{ return gmul(x.g, y.g); }
inline pari gmul2n(pari::tmp x, long n)
{ return gmul2n(x.g, n); }
inline pari gmulsg(long s, pari::tmp y)
{ return gmulsg(s, y.g); }
inline pari gsqr(pari::tmp x)
{ return gsqr(x.g); }
inline pari gsub(pari::tmp x, pari::tmp y)
{ return gsub(x.g, y.g); }
inline pari cgetp(pari::tmp x)
{ return cgetp(x.g); }
inline pari cvtop(pari::tmp x, pari::tmp p, long l)
{ return cvtop(x.g, p.g, l); }
inline pari cvtop2(pari::tmp x, pari::tmp y)
{ return cvtop2(x.g, y.g); }
inline pari gabs(pari::tmp x, long prc = prec)
{ return gabs(x.g, prc); }
inline void gaffect(pari::tmp x, pari::tmp y)
{ gaffect(x.g, y.g); }
inline void gaffsg(long s, pari::tmp x)
{ gaffsg(s, x.g); }
inline int gcmp(pari::tmp x, pari::tmp y)
{ return gcmp(x.g, y.g); }
inline int gcmpsg(long x, pari::tmp y)
{ return gcmpsg(x, y.g); }
inline int gcmp0(pari::tmp x)
{ return gcmp0(x.g); }
inline int gcmp1(pari::tmp x)
{ return gcmp1(x.g); }
inline int gcmp_1(pari::tmp x)
{ return gcmp_1(x.g); }
inline pari gcvtop(pari::tmp x, pari::tmp p, long r)
{ return gcvtop(x.g, p.g, r); }
inline int gequal(pari::tmp x, pari::tmp y)
{ return gequal(x.g, y.g); }
inline int gequalsg(long s, pari::tmp x)
{ return gequalsg(s, x.g); }
inline long gexpo(pari::tmp x)
{ return gexpo(x.g); }
inline long ggval(pari::tmp x, pari::tmp p)
{ return ggval(x.g, p.g); }
inline long glength(pari::tmp x)
{ return glength(x.g); }
inline pari gmax(pari::tmp x, pari::tmp y)
{ return gmax(x.g, y.g); }
inline pari gmaxgs(pari::tmp x, long y)
{ return gmaxgs(x.g, y); }
inline pari gmin(pari::tmp x, pari::tmp y)
{ return gmin(x.g, y.g); }
inline pari gmings(pari::tmp x, long y)
{ return gmings(x.g, y); }
inline pari gneg(pari::tmp x)
{ return gneg(x.g); }
inline pari gneg_i(pari::tmp x)
{ return gneg_i(x.g); }
inline pari greffe(pari::tmp x, long l, long use_stack)
{ return greffe(x.g, l, use_stack); }
inline int gsigne(pari::tmp x)
{ return gsigne(x.g); }
inline pari gtofp(pari::tmp z, long prc = prec)
{ return gtofp(z.g, prc); }
inline pari gtolist(pari::tmp x)
{ return gtolist(x.g); }
inline long gtolong(pari::tmp x)
{ return gtolong(x.g); }
inline int lexcmp(pari::tmp x, pari::tmp y)
{ return lexcmp(x.g, y.g); }
inline pari listconcat(pari::tmp list1, pari::tmp list2)
{ return listconcat(list1.g, list2.g); }
inline pari p_listcreate(long n)
{ return listcreate(n); }
inline pari listinsert(pari::tmp list, pari::tmp object, long index)
{ return listinsert(list.g, object.g, index); }
inline void listkill(pari::tmp list)
{ listkill(list.g); }
inline pari listput(pari::tmp list, pari::tmp object, long index)
{ return listput(list.g, object.g, index); }
inline pari listsort(pari::tmp list, long flag)
{ return listsort(list.g, flag); }
inline pari matsize(pari::tmp x)
{ return matsize(x.g); }
inline pari normalize(pari::tmp x)
{ return normalize(x.g); }
inline pari normalizepol(pari::tmp x)
{ return normalizepol(x.g); }
inline pari normalizepol_approx(pari::tmp x, long lx)
{ return normalizepol_approx(x.g, lx); }
inline pari normalizepol_i(pari::tmp x, long lx)
{ return normalizepol_i(x.g, lx); }
inline pari pureimag(pari::tmp x)
{ return pureimag(x.g); }
inline pari quadtoc(pari::tmp x, long l)
{ return quadtoc(x.g, l); }
inline long sizedigit(pari::tmp x)
{ return sizedigit(x.g); }
inline long u_pvalrem(ulong x, pari::tmp p, ulong * py)
{ return u_pvalrem(x, p.g, py); }
inline pari vecmax(pari::tmp x)
{ return vecmax(x.g); }
inline pari vecmin(pari::tmp x)
{ return vecmin(x.g); }
inline long Z_lval(pari::tmp n, ulong p)
{ return Z_lval(n.g, p); }
inline long z_pval(long n, pari::tmp p)
{ return z_pval(n, p.g); }
inline long Z_pval(pari::tmp n, pari::tmp p)
{ return Z_pval(n.g, p.g); }
inline pari Mod0(pari::tmp x, pari::tmp y, long flag)
{ return Mod0(x.g, y.g, flag); }
inline pari ceil_safe(pari::tmp x)
{ return ceil_safe(x.g); }
inline pari ceilr(pari::tmp x)
{ return ceilr(x.g); }
inline pari centerlift(pari::tmp x)
{ return centerlift(x.g); }
inline pari centerlift0(pari::tmp x, long v)
{ return centerlift0(x.g, v); }
inline pari compo(pari::tmp x, long n)
{ return compo(x.g, n); }
inline pari deg1pol(pari::tmp x1, pari::tmp x0, long v)
{ return deg1pol(x1.g, x0.g, v); }
inline pari deg1pol_i(pari::tmp x1, pari::tmp x0, long v)
{ return deg1pol_i(x1.g, x0.g, v); }
inline long degree(pari::tmp x)
{ return degree(x.g); }
inline pari denom(pari::tmp x)
{ return denom(x.g); }
inline pari deriv(pari::tmp x, long v)
{ return deriv(x.g, v); }
inline pari derivpol(pari::tmp x)
{ return derivpol(x.g); }
inline pari derivser(pari::tmp x)
{ return derivser(x.g); }
inline pari diviiround(pari::tmp x, pari::tmp y)
{ return diviiround(x.g, y.g); }
inline pari divrem(pari::tmp x, pari::tmp y, long v)
{ return divrem(x.g, y.g, v); }
inline pari gand(pari::tmp x, pari::tmp y)
{ return gand(x.g, y.g); }
inline pari gceil(pari::tmp x)
{ return gceil(x.g); }
inline pari gcvtoi(pari::tmp x, long * e)
{ return gcvtoi(x.g, e); }
inline pari gdivent(pari::tmp x, pari::tmp y)
{ return gdivent(x.g, y.g); }
inline pari gdiventgs(pari::tmp x, long y)
{ return gdiventgs(x.g, y); }
inline pari gdiventres(pari::tmp x, pari::tmp y)
{ return gdiventres(x.g, y.g); }
inline pari gdivround(pari::tmp x, pari::tmp y)
{ return gdivround(x.g, y.g); }
inline pari geq(pari::tmp x, pari::tmp y)
{ return geq(x.g, y.g); }
inline pari geval(pari::tmp x)
{ return geval(x.g); }
inline pari gfloor(pari::tmp x)
{ return gfloor(x.g); }
inline pari gfloor2n(pari::tmp x, long s)
{ return gfloor2n(x.g, s); }
inline pari gfrac(pari::tmp x)
{ return gfrac(x.g); }
inline pari gge(pari::tmp x, pari::tmp y)
{ return gge(x.g, y.g); }
inline pari ggprecision(pari::tmp x)
{ return ggprecision(x.g); }
inline pari ggrando(pari::tmp x, long n)
{ return ggrando(x.g, n); }
inline pari ggt(pari::tmp x, pari::tmp y)
{ return ggt(x.g, y.g); }
inline pari gimag(pari::tmp x)
{ return gimag(x.g); }
inline pari ginv(pari::tmp x)
{ return ginv(x.g); }
inline pari gle(pari::tmp x, pari::tmp y)
{ return gle(x.g, y.g); }
inline pari glt(pari::tmp x, pari::tmp y)
{ return glt(x.g, y.g); }
inline pari gmod(pari::tmp x, pari::tmp y)
{ return gmod(x.g, y.g); }
inline pari gmodgs(pari::tmp x, long y)
{ return gmodgs(x.g, y); }
inline pari gmodulo(pari::tmp x, pari::tmp y)
{ return gmodulo(x.g, y.g); }
inline pari gmodulsg(long x, pari::tmp y)
{ return gmodulsg(x, y.g); }
inline pari p_gmodulss(long x, long y)
{ return gmodulss(x, y); }
inline pari gne(pari::tmp x, pari::tmp y)
{ return gne(x.g, y.g); }
inline pari gnot(pari::tmp x)
{ return gnot(x.g); }
inline pari gor(pari::tmp x, pari::tmp y)
{ return gor(x.g, y.g); }
inline pari gpolvar(pari::tmp y)
{ return gpolvar(y.g); }
inline long gprecision(pari::tmp x)
{ return gprecision(x.g); }
inline pari greal(pari::tmp x)
{ return greal(x.g); }
inline pari grndtoi(pari::tmp x, long * e)
{ return grndtoi(x.g, e); }
inline pari ground(pari::tmp x)
{ return ground(x.g); }
inline pari gshift(pari::tmp x, long n)
{ return gshift(x.g, n); }
inline pari gsubst(pari::tmp x, long v, pari::tmp y)
{ return gsubst(x.g, v, y.g); }
inline pari gsubstpol(pari::tmp x, pari::tmp v, pari::tmp y)
{ return gsubstpol(x.g, v.g, y.g); }
inline pari gsubstvec(pari::tmp x, pari::tmp v, pari::tmp y)
{ return gsubstvec(x.g, v.g, y.g); }
inline pari gtocol(pari::tmp x)
{ return gtocol(x.g); }
inline pari gtopoly(pari::tmp x, long v)
{ return gtopoly(x.g, v); }
inline pari gtopolyrev(pari::tmp x, long v)
{ return gtopolyrev(x.g, v); }
inline pari gtoser(pari::tmp x, long v)
{ return gtoser(x.g, v); }
inline pari gtovec(pari::tmp x)
{ return gtovec(x.g); }
inline pari gtovecsmall(pari::tmp x)
{ return gtovecsmall(x.g); }
inline pari gtrunc(pari::tmp x)
{ return gtrunc(x.g); }
inline long gvar(pari::tmp x)
{ return gvar(x.g); }
inline long gvar2(pari::tmp x)
{ return gvar2(x.g); }
inline pari hqfeval(pari::tmp q, pari::tmp x)
{ return hqfeval(q.g, x.g); }
inline pari imag_i(pari::tmp x)
{ return imag_i(x.g); }
inline pari integ(pari::tmp x, long v)
{ return integ(x.g, v); }
inline int iscomplex(pari::tmp x)
{ return iscomplex(x.g); }
inline int isexactzero(pari::tmp g)
{ return isexactzero(g.g); }
inline int isexactzeroscalar(pari::tmp g)
{ return isexactzeroscalar(g.g); }
inline int isinexact(pari::tmp x)
{ return isinexact(x.g); }
inline int isinexactreal(pari::tmp x)
{ return isinexactreal(x.g); }
inline int issmall(pari::tmp n, long * ptk)
{ return issmall(n.g, ptk); }
inline int ismonome(pari::tmp x)
{ return ismonome(x.g); }
inline pari lift(pari::tmp x)
{ return lift(x.g); }
inline pari lift0(pari::tmp x, long v)
{ return lift0(x.g, v); }
inline pari mulmat_real(pari::tmp x, pari::tmp y)
{ return mulmat_real(x.g, y.g); }
inline pari numer(pari::tmp x)
{ return numer(x.g); }
inline long padicprec(pari::tmp x, pari::tmp p)
{ return padicprec(x.g, p.g); }
inline pari polcoeff0(pari::tmp x, long n, long v)
{ return polcoeff0(x.g, n, v); }
inline pari polcoeff_i(pari::tmp x, long n, long v)
{ return polcoeff_i(x.g, n, v); }
inline long poldegree(pari::tmp x, long v)
{ return poldegree(x.g, v); }
inline pari poleval(pari::tmp x, pari::tmp y)
{ return poleval(x.g, y.g); }
inline pari pollead(pari::tmp x, long v)
{ return pollead(x.g, v); }
inline long precision(pari::tmp x)
{ return precision(x.g); }
inline pari precision0(pari::tmp x, long n)
{ return precision0(x.g, n); }
inline pari qf_base_change(pari::tmp q, pari::tmp M, int flag)
{ return qf_base_change(q.g, M.g, flag); }
inline pari qfeval(pari::tmp q, pari::tmp x)
{ return qfeval(q.g, x.g); }
inline pari real_i(pari::tmp x)
{ return real_i(x.g); }
inline pari recip(pari::tmp x)
{ return recip(x.g); }
inline pari roundr(pari::tmp x)
{ return roundr(x.g); }
inline pari scalarpol(pari::tmp x, long v)
{ return scalarpol(x.g, v); }
inline pari scalarser(pari::tmp x, long v, long prc = prec)
{ return scalarser(x.g, v, prc); }
inline pari simplify(pari::tmp x)
{ return simplify(x.g); }
inline pari simplify_i(pari::tmp x)
{ return simplify_i(x.g); }
inline pari tayl(pari::tmp x, long v, long precdl)
{ return tayl(x.g, v, precdl); }
inline pari toser_i(pari::tmp x)
{ return toser_i(x.g); }
inline pari truecoeff(pari::tmp x, long n)
{ return truecoeff(x.g, n); }
inline pari p_u2toi(ulong a, ulong b)
{ return u2toi(a, b); }
inline long group_ident(pari::tmp G, pari::tmp S)
{ return group_ident(G.g, S.g); }
inline long BSW_psp(pari::tmp N)
{ return BSW_psp(N.g); }
inline long millerrabin(pari::tmp n, long k)
{ return millerrabin(n.g, k); }
inline pari nextprime(pari::tmp n)
{ return nextprime(n.g); }
inline pari plisprime(pari::tmp N, long flag)
{ return plisprime(N.g, flag); }
inline pari precprime(pari::tmp n)
{ return precprime(n.g); }
inline pari p_cgetalloc(long t, size_t l)
{ return cgetalloc(t, l); }
inline pari changevar(pari::tmp x, pari::tmp y)
{ return changevar(x.g, y.g); }
inline pari gclone(pari::tmp x)
{ return gclone(x.g); }
inline pari gcopy(pari::tmp x)
{ return gcopy(x.g); }
inline pari gcopy_i(pari::tmp x, long lx)
{ return gcopy_i(x.g, lx); }
inline pari gerepile(pari_sp ltop, pari_sp lbot, pari::tmp q)
{ return gerepile(ltop, lbot, q.g); }
inline void gerepilecoeffs(pari_sp av, pari::tmp x, int n)
{ gerepilecoeffs(av, x.g, n); }
inline pari gerepilecopy(pari_sp av, pari::tmp x)
{ return gerepilecopy(av, x.g); }
inline pari gerepileupto(pari_sp av, pari::tmp q)
{ return gerepileupto(av, q.g); }
inline pari gerepileuptoint(pari_sp av, pari::tmp q)
{ return gerepileuptoint(av, q.g); }
inline pari gerepileuptoleaf(pari_sp av, pari::tmp q)
{ return gerepileuptoleaf(av, q.g); }
inline void gunclone(pari::tmp x)
{ gunclone(x.g); }
inline void killbloc(pari::tmp x)
{ killbloc(x.g); }
inline pari p_newbloc(long n)
{ return newbloc(n); }
inline pari reorder(pari::tmp x)
{ return reorder(x.g); }
inline pari shallowcopy(pari::tmp x)
{ return shallowcopy(x.g); }
inline long taille(pari::tmp x)
{ return taille(x.g); }
inline long taille2(pari::tmp x)
{ return taille2(x.g); }
inline pari intcirc0(entree * ep, pari::tmp a, pari::tmp R, char * ch, pari::tmp tab, long prc = prec)
{ return intcirc0(ep, a.g, R.g, ch, tab.g, prc); }
inline pari intfourcos0(entree * ep, pari::tmp a, pari::tmp b, pari::tmp x, char * ch, pari::tmp tab, long prc = prec)
{ return intfourcos0(ep, a.g, b.g, x.g, ch, tab.g, prc); }
inline pari intfourexp0(entree * ep, pari::tmp a, pari::tmp b, pari::tmp x, char * ch, pari::tmp tab, long prc = prec)
{ return intfourexp0(ep, a.g, b.g, x.g, ch, tab.g, prc); }
inline pari intfoursin0(entree * ep, pari::tmp a, pari::tmp b, pari::tmp x, char * ch, pari::tmp tab, long prc = prec)
{ return intfoursin0(ep, a.g, b.g, x.g, ch, tab.g, prc); }
inline pari intfuncinit0(entree * ep, pari::tmp a, pari::tmp b, char * ch, long flag, long m, long prc = prec)
{ return intfuncinit0(ep, a.g, b.g, ch, flag, m, prc); }
inline pari intlaplaceinv0(entree * ep, pari::tmp sig, pari::tmp x, char * ch, pari::tmp tab, long prc = prec)
{ return intlaplaceinv0(ep, sig.g, x.g, ch, tab.g, prc); }
inline pari intmellininv0(entree * ep, pari::tmp sig, pari::tmp x, char * ch, pari::tmp tab, long prc = prec)
{ return intmellininv0(ep, sig.g, x.g, ch, tab.g, prc); }
inline pari intmellininvshort(pari::tmp sig, pari::tmp x, pari::tmp tab, long prc = prec)
{ return intmellininvshort(sig.g, x.g, tab.g, prc); }
inline pari intnum0(entree * ep, pari::tmp a, pari::tmp b, char * ch, pari::tmp tab, long prc = prec)
{ return intnum0(ep, a.g, b.g, ch, tab.g, prc); }
inline pari intnuminit(pari::tmp a, pari::tmp b, long m, long prc = prec)
{ return intnuminit(a.g, b.g, m, prc); }
inline pari intnuminit0(pari::tmp a, pari::tmp b, pari::tmp tab, long prc = prec)
{ return intnuminit0(a.g, b.g, tab.g, prc); }
inline pari intnuminitgen0(entree * ep, pari::tmp a, pari::tmp b, char * ch, long m, long flag, long prc = prec)
{ return intnuminitgen0(ep, a.g, b.g, ch, m, flag, prc); }
inline pari intnumromb0(entree * ep, pari::tmp a, pari::tmp b, char * ch, long flag, long prc = prec)
{ return intnumromb0(ep, a.g, b.g, ch, flag, prc); }
inline pari sumnum0(entree * ep, pari::tmp a, pari::tmp sig, char * ch, pari::tmp tab, long flag, long prc = prec)
{ return sumnum0(ep, a.g, sig.g, ch, tab.g, flag, prc); }
inline pari sumnumalt0(entree * ep, pari::tmp a, pari::tmp sig, char * ch, pari::tmp tab, long flag, long prc = prec)
{ return sumnumalt0(ep, a.g, sig.g, ch, tab.g, flag, prc); }
inline pari sumnuminit(pari::tmp sig, long m, long sgn, long prc = prec)
{ return sumnuminit(sig.g, m, sgn, prc); }
inline pari sumnuminit0(pari::tmp a, pari::tmp tab, long sgn, long prc = prec)
{ return sumnuminit0(a.g, tab.g, sgn, prc); }
inline pari rnfkummer(pari::tmp bnr, pari::tmp subgroup, long all, long prc = prec)
{ return rnfkummer(bnr.g, subgroup.g, all, prc); }
inline pari member_a1(pari::tmp x)
{ return member_a1(x.g); }
inline pari member_a2(pari::tmp x)
{ return member_a2(x.g); }
inline pari member_a3(pari::tmp x)
{ return member_a3(x.g); }
inline pari member_a4(pari::tmp x)
{ return member_a4(x.g); }
inline pari member_a6(pari::tmp x)
{ return member_a6(x.g); }
inline pari member_area(pari::tmp x)
{ return member_area(x.g); }
inline pari member_b2(pari::tmp x)
{ return member_b2(x.g); }
inline pari member_b4(pari::tmp x)
{ return member_b4(x.g); }
inline pari member_b6(pari::tmp x)
{ return member_b6(x.g); }
inline pari member_b8(pari::tmp x)
{ return member_b8(x.g); }
inline pari member_bid(pari::tmp x)
{ return member_bid(x.g); }
inline pari member_bnf(pari::tmp x)
{ return member_bnf(x.g); }
inline pari member_c4(pari::tmp x)
{ return member_c4(x.g); }
inline pari member_c6(pari::tmp x)
{ return member_c6(x.g); }
inline pari member_clgp(pari::tmp x)
{ return member_clgp(x.g); }
inline pari member_codiff(pari::tmp x)
{ return member_codiff(x.g); }
inline pari member_cyc(pari::tmp clg)
{ return member_cyc(clg.g); }
inline pari member_diff(pari::tmp x)
{ return member_diff(x.g); }
inline pari member_disc(pari::tmp x)
{ return member_disc(x.g); }
inline pari member_e(pari::tmp x)
{ return member_e(x.g); }
inline pari member_eta(pari::tmp x)
{ return member_eta(x.g); }
inline pari member_f(pari::tmp x)
{ return member_f(x.g); }
inline pari member_fu(pari::tmp x)
{ return member_fu(x.g); }
inline pari member_futu(pari::tmp x)
{ return member_futu(x.g); }
inline pari member_gen(pari::tmp x)
{ return member_gen(x.g); }
inline pari member_group(pari::tmp x)
{ return member_group(x.g); }
inline pari member_index(pari::tmp x)
{ return member_index(x.g); }
inline pari member_j(pari::tmp x)
{ return member_j(x.g); }
inline pari member_mod(pari::tmp x)
{ return member_mod(x.g); }
inline pari member_nf(pari::tmp x)
{ return member_nf(x.g); }
inline pari member_no(pari::tmp clg)
{ return member_no(clg.g); }
inline pari member_omega(pari::tmp x)
{ return member_omega(x.g); }
inline pari member_orders(pari::tmp x)
{ return member_orders(x.g); }
inline pari member_p(pari::tmp x)
{ return member_p(x.g); }
inline pari member_pol(pari::tmp x)
{ return member_pol(x.g); }
inline pari member_reg(pari::tmp x)
{ return member_reg(x.g); }
inline pari member_r1(pari::tmp x)
{ return member_r1(x.g); }
inline pari member_r2(pari::tmp x)
{ return member_r2(x.g); }
inline pari member_roots(pari::tmp x)
{ return member_roots(x.g); }
inline pari member_sign(pari::tmp x)
{ return member_sign(x.g); }
inline pari member_t2(pari::tmp x)
{ return member_t2(x.g); }
inline pari member_tate(pari::tmp x)
{ return member_tate(x.g); }
inline pari member_tufu(pari::tmp x)
{ return member_tufu(x.g); }
inline pari member_tu(pari::tmp x)
{ return member_tu(x.g); }
inline pari member_w(pari::tmp x)
{ return member_w(x.g); }
inline pari member_zk(pari::tmp x)
{ return member_zk(x.g); }
inline pari member_zkst(pari::tmp bid)
{ return member_zkst(bid.g); }
inline int absi_cmp(pari::tmp x, pari::tmp y)
{ return absi_cmp(x.g, y.g); }
inline int absi_equal(pari::tmp x, pari::tmp y)
{ return absi_equal(x.g, y.g); }
inline int absr_cmp(pari::tmp x, pari::tmp y)
{ return absr_cmp(x.g, y.g); }
inline pari addii_sign(pari::tmp x, long sx, pari::tmp y, long sy)
{ return addii_sign(x.g, sx, y.g, sy); }
inline pari addir_sign(pari::tmp x, long sx, pari::tmp y, long sy)
{ return addir_sign(x.g, sx, y.g, sy); }
inline pari addrr_sign(pari::tmp x, long sx, pari::tmp y, long sy)
{ return addrr_sign(x.g, sx, y.g, sy); }
inline pari addsi_sign(long x, pari::tmp y, long sy)
{ return addsi_sign(x, y.g, sy); }
inline pari addsr(long x, pari::tmp y)
{ return addsr(x, y.g); }
inline pari p_addss(long x, long y)
{ return addss(x, y); }
inline void affir(pari::tmp x, pari::tmp y)
{ affir(x.g, y.g); }
inline void affrr(pari::tmp x, pari::tmp y)
{ affrr(x.g, y.g); }
inline void cgiv(pari::tmp x)
{ cgiv(x.g); }
inline int cmpii(pari::tmp x, pari::tmp y)
{ return cmpii(x.g, y.g); }
inline int cmprr(pari::tmp x, pari::tmp y)
{ return cmprr(x.g, y.g); }
inline int cmpsi(long x, pari::tmp y)
{ return cmpsi(x, y.g); }
inline int cmpui(ulong x, pari::tmp y)
{ return cmpui(x, y.g); }
inline pari p_dbltor(double x)
{ return dbltor(x); }
inline pari diviiexact(pari::tmp x, pari::tmp y)
{ return diviiexact(x.g, y.g); }
inline pari diviuexact(pari::tmp x, ulong y)
{ return diviuexact(x.g, y); }
inline pari divir(pari::tmp x, pari::tmp y)
{ return divir(x.g, y.g); }
inline pari divis(pari::tmp y, long x)
{ return divis(y.g, x); }
inline pari divis_rem(pari::tmp x, long y, long * rem)
{ return divis_rem(x.g, y, rem); }
inline pari diviu_rem(pari::tmp y, ulong x, ulong * rem)
{ return diviu_rem(y.g, x, rem); }
inline pari divri(pari::tmp x, pari::tmp y)
{ return divri(x.g, y.g); }
inline pari divrr(pari::tmp x, pari::tmp y)
{ return divrr(x.g, y.g); }
inline pari divrs(pari::tmp x, long y)
{ return divrs(x.g, y); }
inline pari divsi(long x, pari::tmp y)
{ return divsi(x, y.g); }
inline pari divsr(long x, pari::tmp y)
{ return divsr(x, y.g); }
inline int equalii(pari::tmp x, pari::tmp y)
{ return equalii(x.g, y.g); }
inline int equalsi(long x, pari::tmp y)
{ return equalsi(x, y.g); }
inline int equalui(ulong x, pari::tmp y)
{ return equalui(x, y.g); }
inline pari floorr(pari::tmp x)
{ return floorr(x.g); }
inline pari gcdii(pari::tmp x, pari::tmp y)
{ return gcdii(x.g, y.g); }
inline pari int_normalize(pari::tmp x, long known_zero_words)
{ return int_normalize(x.g, known_zero_words); }
inline pari p_int2n(long n)
{ return int2n(n); }
inline pari p_int2u(ulong n)
{ return int2u(n); }
inline pari ishiftr(pari::tmp x, long n)
{ return ishiftr(x.g, n); }
inline pari modii(pari::tmp x, pari::tmp y)
{ return modii(x.g, y.g); }
inline void modiiz(pari::tmp x, pari::tmp y, pari::tmp z)
{ modiiz(x.g, y.g, z.g); }
inline void mpdivz(pari::tmp x, pari::tmp y, pari::tmp z)
{ mpdivz(x.g, y.g, z.g); }
inline pari mulii(pari::tmp x, pari::tmp y)
{ return mulii(x.g, y.g); }
inline pari mulir(pari::tmp x, pari::tmp y)
{ return mulir(x.g, y.g); }
inline pari mulrr(pari::tmp x, pari::tmp y)
{ return mulrr(x.g, y.g); }
inline pari mulsi(long x, pari::tmp y)
{ return mulsi(x, y.g); }
inline pari mulsr(long x, pari::tmp y)
{ return mulsr(x, y.g); }
inline pari p_mulss(long x, long y)
{ return mulss(x, y); }
inline pari mului(ulong x, pari::tmp y)
{ return mului(x, y.g); }
inline pari mulur(ulong x, pari::tmp y)
{ return mulur(x, y.g); }
inline pari p_muluu(ulong x, ulong y)
{ return muluu(x, y); }
inline pari randomi(pari::tmp x)
{ return randomi(x.g); }
inline pari resmod2n(pari::tmp x, long n)
{ return resmod2n(x.g, n); }
inline double rtodbl(pari::tmp x)
{ return rtodbl(x.g); }
inline pari shifti(pari::tmp x, long n)
{ return shifti(x.g, n); }
inline pari sqri(pari::tmp x)
{ return sqri(x.g); }
inline pari subsr(long x, pari::tmp y)
{ return subsr(x, y.g); }
inline pari truncr(pari::tmp x)
{ return truncr(x.g); }
inline ulong umodiu(pari::tmp y, ulong x)
{ return umodiu(y.g, x); }
inline pari nffactor(pari::tmp nf, pari::tmp x)
{ return nffactor(nf.g, x.g); }
inline pari nffactormod(pari::tmp nf, pari::tmp pol, pari::tmp pr)
{ return nffactormod(nf.g, pol.g, pr.g); }
inline int nfisgalois(pari::tmp nf, pari::tmp x)
{ return nfisgalois(nf.g, x.g); }
inline pari nfroots(pari::tmp nf, pari::tmp pol)
{ return nfroots(nf.g, pol.g); }
inline pari rnfcharpoly(pari::tmp nf, pari::tmp T, pari::tmp alpha, long v)
{ return rnfcharpoly(nf.g, T.g, alpha.g, v); }
inline pari unifpol(pari::tmp nf, pari::tmp pol, long flag)
{ return unifpol(nf.g, pol.g, flag); }
inline pari numbpart(pari::tmp x)
{ return numbpart(x.g); }
inline pari abelian_group(pari::tmp G)
{ return abelian_group(G.g); }
inline pari p_bitvec_alloc(long n)
{ return bitvec_alloc(n); }
inline void bitvec_clear(pari::tmp bitvec, long b)
{ bitvec_clear(bitvec.g, b); }
inline void bitvec_set(pari::tmp bitvec, long b)
{ bitvec_set(bitvec.g, b); }
inline pari bitvec_shorten(pari::tmp bitvec, long n)
{ return bitvec_shorten(bitvec.g, n); }
inline long bitvec_test(pari::tmp bitvec, long b)
{ return bitvec_test(bitvec.g, b); }
inline pari const_col(long n, pari::tmp x)
{ return const_col(n, x.g); }
inline pari const_vec(long n, pari::tmp x)
{ return const_vec(n, x.g); }
inline pari p_const_vecsmall(long n, long c)
{ return const_vecsmall(n, c); }
inline pari cyclicgroup(pari::tmp g, long s)
{ return cyclicgroup(g.g, s); }
inline pari p_cyclicperm(long l, long d)
{ return cyclicperm(l, d); }
inline pari cyc_pow(pari::tmp cyc, long exp)
{ return cyc_pow(cyc.g, exp); }
inline pari cyc_pow_perm(pari::tmp cyc, long exp)
{ return cyc_pow_perm(cyc.g, exp); }
inline pari dicyclicgroup(pari::tmp g1, pari::tmp g2, long s1, long s2)
{ return dicyclicgroup(g1.g, g2.g, s1, s2); }
inline pari group_abelianHNF(pari::tmp G, pari::tmp L)
{ return group_abelianHNF(G.g, L.g); }
inline pari group_abelianSNF(pari::tmp G, pari::tmp L)
{ return group_abelianSNF(G.g, L.g); }
inline long group_domain(pari::tmp G)
{ return group_domain(G.g); }
inline pari group_elts(pari::tmp G, long n)
{ return group_elts(G.g, n); }
inline pari group_export(pari::tmp G, long format)
{ return group_export(G.g, format); }
inline long group_isA4S4(pari::tmp G)
{ return group_isA4S4(G.g); }
inline long group_isabelian(pari::tmp G)
{ return group_isabelian(G.g); }
inline pari group_leftcoset(pari::tmp G, pari::tmp g)
{ return group_leftcoset(G.g, g.g); }
inline long group_order(pari::tmp G)
{ return group_order(G.g); }
inline long group_perm_normalize(pari::tmp N, pari::tmp g)
{ return group_perm_normalize(N.g, g.g); }
inline pari group_quotient(pari::tmp G, pari::tmp H)
{ return group_quotient(G.g, H.g); }
inline pari group_rightcoset(pari::tmp G, pari::tmp g)
{ return group_rightcoset(G.g, g.g); }
inline pari group_subgroups(pari::tmp G)
{ return group_subgroups(G.g); }
inline pari groupelts_center(pari::tmp S)
{ return groupelts_center(S.g); }
inline pari groupelts_abelian_group(pari::tmp S)
{ return groupelts_abelian_group(S.g); }
inline int perm_commute(pari::tmp p, pari::tmp q)
{ return perm_commute(p.g, q.g); }
inline pari perm_cycles(pari::tmp v)
{ return perm_cycles(v.g); }
inline pari p_perm_identity(long l)
{ return perm_identity(l); }
inline pari perm_inv(pari::tmp x)
{ return perm_inv(x.g); }
inline pari perm_mul(pari::tmp s, pari::tmp t)
{ return perm_mul(s.g, t.g); }
inline long perm_order(pari::tmp perm)
{ return perm_order(perm.g); }
inline pari perm_pow(pari::tmp perm, long exp)
{ return perm_pow(perm.g, exp); }
inline pari quotient_group(pari::tmp C, pari::tmp G)
{ return quotient_group(C.g, G.g); }
inline pari quotient_perm(pari::tmp C, pari::tmp p)
{ return quotient_perm(C.g, p.g); }
inline pari vec_to_vecsmall(pari::tmp z)
{ return vec_to_vecsmall(z.g); }
inline pari vecperm_orbits(pari::tmp v, long n)
{ return vecperm_orbits(v.g, n); }
inline int vec_is1to1(pari::tmp v)
{ return vec_is1to1(v.g); }
inline int vec_isconst(pari::tmp v)
{ return vec_isconst(v.g); }
inline pari vec_lengthen(pari::tmp v, long n)
{ return vec_lengthen(v.g, n); }
inline pari vec_shorten(pari::tmp v, long n)
{ return vec_shorten(v.g, n); }
inline pari vecsmall_append(pari::tmp V, long s)
{ return vecsmall_append(V.g, s); }
inline long vecsmall_coincidence(pari::tmp u, pari::tmp v)
{ return vecsmall_coincidence(u.g, v.g); }
inline pari vecsmall_concat(pari::tmp u, pari::tmp v)
{ return vecsmall_concat(u.g, v.g); }
inline pari vecsmall_copy(pari::tmp x)
{ return vecsmall_copy(x.g); }
inline pari vecsmall_indexsort(pari::tmp V)
{ return vecsmall_indexsort(V.g); }
inline long vecsmall_isin(pari::tmp v, long x)
{ return vecsmall_isin(v.g, x); }
inline pari vecsmall_lengthen(pari::tmp v, long n)
{ return vecsmall_lengthen(v.g, n); }
inline int vecsmall_lexcmp(pari::tmp x, pari::tmp y)
{ return vecsmall_lexcmp(x.g, y.g); }
inline long vecsmall_pack(pari::tmp V, long base, long mod)
{ return vecsmall_pack(V.g, base, mod); }
inline int vecsmall_prefixcmp(pari::tmp x, pari::tmp y)
{ return vecsmall_prefixcmp(x.g, y.g); }
inline pari vecsmall_prepend(pari::tmp V, long s)
{ return vecsmall_prepend(V.g, s); }
inline pari vecsmall_shorten(pari::tmp v, long n)
{ return vecsmall_shorten(v.g, n); }
inline void vecsmall_sort(pari::tmp V)
{ vecsmall_sort(V.g); }
inline pari vecsmall_to_col(pari::tmp z)
{ return vecsmall_to_col(z.g); }
inline pari vecsmall_to_vec(pari::tmp z)
{ return vecsmall_to_vec(z.g); }
inline pari vecsmall_uniq(pari::tmp V)
{ return vecsmall_uniq(V.g); }
inline pari vecvecsmall_indexsort(pari::tmp x)
{ return vecvecsmall_indexsort(x.g); }
inline pari vecvecsmall_sort(pari::tmp x)
{ return vecvecsmall_sort(x.g); }
inline long vecvecsmall_search(pari::tmp x, pari::tmp y, long flag)
{ return vecvecsmall_search(x.g, y.g, flag); }
inline long Flx_nbfact(pari::tmp z, ulong p)
{ return Flx_nbfact(z.g, p); }
inline long Flx_nbroots(pari::tmp f, ulong p)
{ return Flx_nbroots(f.g, p); }
inline pari FpX_degfact(pari::tmp f, pari::tmp p)
{ return FpX_degfact(f.g, p.g); }
inline long FpX_is_irred(pari::tmp f, pari::tmp p)
{ return FpX_is_irred(f.g, p.g); }
inline long FpX_is_squarefree(pari::tmp f, pari::tmp p)
{ return FpX_is_squarefree(f.g, p.g); }
inline long FpX_is_totally_split(pari::tmp f, pari::tmp p)
{ return FpX_is_totally_split(f.g, p.g); }
inline pari FpX_factor(pari::tmp f, pari::tmp p)
{ return FpX_factor(f.g, p.g); }
inline long FpX_nbfact(pari::tmp f, pari::tmp p)
{ return FpX_nbfact(f.g, p.g); }
inline long FpX_nbroots(pari::tmp f, pari::tmp p)
{ return FpX_nbroots(f.g, p.g); }
inline pari FqX_factor(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FqX_factor(x.g, T.g, p.g); }
inline pari FqX_gcd(pari::tmp P, pari::tmp Q, pari::tmp T, pari::tmp p)
{ return FqX_gcd(P.g, Q.g, T.g, p.g); }
inline long FqX_is_squarefree(pari::tmp P, pari::tmp T, pari::tmp p)
{ return FqX_is_squarefree(P.g, T.g, p.g); }
inline long FqX_nbfact(pari::tmp u, pari::tmp T, pari::tmp p)
{ return FqX_nbfact(u.g, T.g, p.g); }
inline long FqX_nbroots(pari::tmp f, pari::tmp T, pari::tmp p)
{ return FqX_nbroots(f.g, T.g, p.g); }
inline pari FpX_rand(long d, long v, pari::tmp p)
{ return FpX_rand(d, v, p.g); }
inline pari FpX_roots(pari::tmp f, pari::tmp p)
{ return FpX_roots(f.g, p.g); }
inline int cmp_pol(pari::tmp x, pari::tmp y)
{ return cmp_pol(x.g, y.g); }
inline pari factcantor(pari::tmp x, pari::tmp p)
{ return factcantor(x.g, p.g); }
inline pari factmod(pari::tmp f, pari::tmp p)
{ return factmod(f.g, p.g); }
inline pari factorff(pari::tmp f, pari::tmp p, pari::tmp a)
{ return factorff(f.g, p.g, a.g); }
inline pari factormod0(pari::tmp f, pari::tmp p, long flag)
{ return factormod0(f.g, p.g, flag); }
inline pari factorpadic0(pari::tmp f, pari::tmp p, long r, long flag)
{ return factorpadic0(f.g, p.g, r, flag); }
inline pari factorpadic2(pari::tmp x, pari::tmp p, long r)
{ return factorpadic2(x.g, p.g, r); }
inline pari factorpadic4(pari::tmp x, pari::tmp p, long r)
{ return factorpadic4(x.g, p.g, r); }
inline int gdvd(pari::tmp x, pari::tmp y)
{ return gdvd(x.g, y.g); }
inline pari padicappr(pari::tmp f, pari::tmp a)
{ return padicappr(f.g, a.g); }
inline pari padicsqrtnlift(pari::tmp a, pari::tmp n, pari::tmp S, pari::tmp p, long e)
{ return padicsqrtnlift(a.g, n.g, S.g, p.g, e); }
inline pari polfnf(pari::tmp a, pari::tmp t)
{ return polfnf(a.g, t.g); }
inline pari rootmod(pari::tmp f, pari::tmp p)
{ return rootmod(f.g, p.g); }
inline pari rootmod0(pari::tmp f, pari::tmp p, long flag)
{ return rootmod0(f.g, p.g, flag); }
inline pari rootmod2(pari::tmp f, pari::tmp p)
{ return rootmod2(f.g, p.g); }
inline pari rootpadic(pari::tmp f, pari::tmp p, long r)
{ return rootpadic(f.g, p.g, r); }
inline pari rootpadicfast(pari::tmp f, pari::tmp p, long e)
{ return rootpadicfast(f.g, p.g, e); }
inline pari ZX_deriv(pari::tmp x)
{ return ZX_deriv(x.g); }
inline pari FpX_deriv(pari::tmp f, pari::tmp p)
{ return FpX_deriv(f.g, p.g); }
inline pari FqX_deriv(pari::tmp f, pari::tmp T, pari::tmp p)
{ return FqX_deriv(f.g, T.g, p.g); }
inline pari ZpX_liftroot(pari::tmp f, pari::tmp a, pari::tmp p, long e)
{ return ZpX_liftroot(f.g, a.g, p.g, e); }
inline pari ZpXQX_liftroot(pari::tmp f, pari::tmp a, pari::tmp T, pari::tmp p, long e)
{ return ZpXQX_liftroot(f.g, a.g, T.g, p.g, e); }
inline pari ZpX_liftroots(pari::tmp f, pari::tmp S, pari::tmp q, long e)
{ return ZpX_liftroots(f.g, S.g, q.g, e); }
inline pari roots2(pari::tmp pol, long PREC)
{ return roots2(pol.g, PREC); }
inline pari rootsold(pari::tmp x, long l)
{ return rootsold(x.g, l); }
inline pari simplefactmod(pari::tmp f, pari::tmp p)
{ return simplefactmod(f.g, p.g); }
inline pari p_Newton_exponents(long e)
{ return Newton_exponents(e); }
inline pari Q_content(pari::tmp x)
{ return Q_content(x.g); }
inline pari Q_denom(pari::tmp x)
{ return Q_denom(x.g); }
inline pari Q_div_to_int(pari::tmp x, pari::tmp c)
{ return Q_div_to_int(x.g, c.g); }
inline pari Q_muli_to_int(pari::tmp x, pari::tmp d)
{ return Q_muli_to_int(x.g, d.g); }
inline pari Q_primpart(pari::tmp x)
{ return Q_primpart(x.g); }
inline pari centermod(pari::tmp x, pari::tmp p)
{ return centermod(x.g, p.g); }
inline pari centermod_i(pari::tmp x, pari::tmp p, pari::tmp ps2)
{ return centermod_i(x.g, p.g, ps2.g); }
inline pari centermodii(pari::tmp x, pari::tmp p, pari::tmp po2)
{ return centermodii(x.g, p.g, po2.g); }
inline pari combine_factors(pari::tmp target, pari::tmp famod, pari::tmp p, long klim, long hint)
{ return combine_factors(target.g, famod.g, p.g, klim, hint); }
inline pari concat_factor(pari::tmp f, pari::tmp g)
{ return concat_factor(f.g, g.g); }
inline pari content(pari::tmp x)
{ return content(x.g); }
inline pari deg1_from_roots(pari::tmp L, long v)
{ return deg1_from_roots(L.g, v); }
inline pari discsr(pari::tmp x)
{ return discsr(x.g); }
inline pari factor(pari::tmp x)
{ return factor(x.g); }
inline pari factor0(pari::tmp x, long flag)
{ return factor0(x.g, flag); }
inline pari factorback(pari::tmp fa, pari::tmp nf)
{ return factorback(fa.g, nf.g); }
inline pari factorback0(pari::tmp fa, pari::tmp e, pari::tmp nf)
{ return factorback0(fa.g, e.g, nf.g); }
inline pari factorbackelt(pari::tmp fa, pari::tmp e, pari::tmp nf)
{ return factorbackelt(fa.g, e.g, nf.g); }
inline pari factpol(pari::tmp x, long hint)
{ return factpol(x.g, hint); }
inline pari gcd0(pari::tmp x, pari::tmp y, long flag)
{ return gcd0(x.g, y.g, flag); }
inline pari gdeflate(pari::tmp x, long v, long d)
{ return gdeflate(x.g, v, d); }
inline pari gdivexact(pari::tmp x, pari::tmp y)
{ return gdivexact(x.g, y.g); }
inline pari ggcd(pari::tmp x, pari::tmp y)
{ return ggcd(x.g, y.g); }
inline pari ginvmod(pari::tmp x, pari::tmp y)
{ return ginvmod(x.g, y.g); }
inline pari gisirreducible(pari::tmp x)
{ return gisirreducible(x.g); }
inline pari glcm(pari::tmp x, pari::tmp y)
{ return glcm(x.g, y.g); }
inline pari glcm0(pari::tmp x, pari::tmp y)
{ return glcm0(x.g, y.g); }
inline pari hensel_lift_fact(pari::tmp pol, pari::tmp Q, pari::tmp T, pari::tmp p, pari::tmp pe, long e)
{ return hensel_lift_fact(pol.g, Q.g, T.g, p.g, pe.g, e); }
inline pari newtonpoly(pari::tmp x, pari::tmp p)
{ return newtonpoly(x.g, p.g); }
inline pari nfgcd(pari::tmp P, pari::tmp Q, pari::tmp nf, pari::tmp den)
{ return nfgcd(P.g, Q.g, nf.g, den.g); }
inline pari nfrootsQ(pari::tmp x)
{ return nfrootsQ(x.g); }
inline pari poldeflate(pari::tmp x0, long * m)
{ return poldeflate(x0.g, m); }
inline pari poldeflate_i(pari::tmp x0, long d)
{ return poldeflate_i(x0.g, d); }
inline pari poldisc0(pari::tmp x, long v)
{ return poldisc0(x.g, v); }
inline pari polhensellift(pari::tmp pol, pari::tmp fct, pari::tmp p, long exp)
{ return polhensellift(pol.g, fct.g, p.g, exp); }
inline pari polinflate(pari::tmp x0, long d)
{ return polinflate(x0.g, d); }
inline pari polresultant0(pari::tmp x, pari::tmp y, long v, long flag)
{ return polresultant0(x.g, y.g, v, flag); }
inline pari polsym(pari::tmp x, long n)
{ return polsym(x.g, n); }
inline pari primpart(pari::tmp x)
{ return primpart(x.g); }
inline pari pseudorem(pari::tmp x, pari::tmp y)
{ return pseudorem(x.g, y.g); }
inline pari reduceddiscsmith(pari::tmp pol)
{ return reduceddiscsmith(pol.g); }
inline pari resultant2(pari::tmp x, pari::tmp y)
{ return resultant2(x.g, y.g); }
inline pari resultantducos(pari::tmp x, pari::tmp y)
{ return resultantducos(x.g, y.g); }
inline pari roots_from_deg1(pari::tmp x)
{ return roots_from_deg1(x.g); }
inline pari srgcd(pari::tmp x, pari::tmp y)
{ return srgcd(x.g, y.g); }
inline long sturmpart(pari::tmp x, pari::tmp a, pari::tmp b)
{ return sturmpart(x.g, a.g, b.g); }
inline pari sylvestermatrix(pari::tmp x, pari::tmp y)
{ return sylvestermatrix(x.g, y.g); }
inline pari vecbezout(pari::tmp x, pari::tmp y)
{ return vecbezout(x.g, y.g); }
inline pari vecbezoutres(pari::tmp x, pari::tmp y)
{ return vecbezoutres(x.g, y.g); }
inline pari FpC_red(pari::tmp z, pari::tmp p)
{ return FpC_red(z.g, p.g); }
inline pari FpC_to_mod(pari::tmp z, pari::tmp p)
{ return FpC_to_mod(z.g, p.g); }
inline pari FpM_red(pari::tmp z, pari::tmp p)
{ return FpM_red(z.g, p.g); }
inline pari FpM_to_mod(pari::tmp z, pari::tmp p)
{ return FpM_to_mod(z.g, p.g); }
inline pari FpV_polint(pari::tmp xa, pari::tmp ya, pari::tmp p)
{ return FpV_polint(xa.g, ya.g, p.g); }
inline pari FpV_red(pari::tmp z, pari::tmp p)
{ return FpV_red(z.g, p.g); }
inline pari FpV_roots_to_pol(pari::tmp V, pari::tmp p, long v)
{ return FpV_roots_to_pol(V.g, p.g, v); }
inline pari FpV_to_mod(pari::tmp z, pari::tmp p)
{ return FpV_to_mod(z.g, p.g); }
inline pari FpX_FpXQ_compo(pari::tmp f, pari::tmp x, pari::tmp T, pari::tmp p)
{ return FpX_FpXQ_compo(f.g, x.g, T.g, p.g); }
inline pari FpX_FpXQV_compo(pari::tmp f, pari::tmp x, pari::tmp T, pari::tmp p)
{ return FpX_FpXQV_compo(f.g, x.g, T.g, p.g); }
inline pari FpX_Fp_add(pari::tmp y, pari::tmp x, pari::tmp p)
{ return FpX_Fp_add(y.g, x.g, p.g); }
inline pari FpX_Fp_mul(pari::tmp y, pari::tmp x, pari::tmp p)
{ return FpX_Fp_mul(y.g, x.g, p.g); }
inline pari FpX_add(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpX_add(x.g, y.g, p.g); }
inline pari FpX_center(pari::tmp T, pari::tmp mod)
{ return FpX_center(T.g, mod.g); }
inline pari FpX_chinese_coprime(pari::tmp x, pari::tmp y, pari::tmp Tx, pari::tmp Ty, pari::tmp Tz, pari::tmp p)
{ return FpX_chinese_coprime(x.g, y.g, Tx.g, Ty.g, Tz.g, p.g); }
inline pari FpX_eval(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpX_eval(x.g, y.g, p.g); }
inline pari FpX_factorff_irred(pari::tmp P, pari::tmp Q, pari::tmp l)
{ return FpX_factorff_irred(P.g, Q.g, l.g); }
inline pari FpX_ffisom(pari::tmp P, pari::tmp Q, pari::tmp l)
{ return FpX_ffisom(P.g, Q.g, l.g); }
inline pari FpX_gcd(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpX_gcd(x.g, y.g, p.g); }
inline pari FpX_mul(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpX_mul(x.g, y.g, p.g); }
inline pari FpX_neg(pari::tmp x, pari::tmp p)
{ return FpX_neg(x.g, p.g); }
inline pari FpX_normalize(pari::tmp z, pari::tmp p)
{ return FpX_normalize(z.g, p.g); }
inline pari FpX_red(pari::tmp z, pari::tmp p)
{ return FpX_red(z.g, p.g); }
inline pari FpX_resultant(pari::tmp a, pari::tmp b, pari::tmp p)
{ return FpX_resultant(a.g, b.g, p.g); }
inline pari FpX_sqr(pari::tmp x, pari::tmp p)
{ return FpX_sqr(x.g, p.g); }
inline pari FpX_sub(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpX_sub(x.g, y.g, p.g); }
inline pari FpX_to_mod(pari::tmp z, pari::tmp p)
{ return FpX_to_mod(z.g, p.g); }
inline pari FpXQ_charpoly(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FpXQ_charpoly(x.g, T.g, p.g); }
inline pari FpXQ_div(pari::tmp x, pari::tmp y, pari::tmp T, pari::tmp p)
{ return FpXQ_div(x.g, y.g, T.g, p.g); }
inline pari FpXQ_ffisom_inv(pari::tmp S, pari::tmp Tp, pari::tmp p)
{ return FpXQ_ffisom_inv(S.g, Tp.g, p.g); }
inline pari FpXQ_inv(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FpXQ_inv(x.g, T.g, p.g); }
inline pari FpXQ_invsafe(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FpXQ_invsafe(x.g, T.g, p.g); }
inline pari FpXQ_matrix_pow(pari::tmp y, long n, long m, pari::tmp P, pari::tmp l)
{ return FpXQ_matrix_pow(y.g, n, m, P.g, l.g); }
inline pari FpXQ_minpoly(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FpXQ_minpoly(x.g, T.g, p.g); }
inline pari FpXQ_mul(pari::tmp y, pari::tmp x, pari::tmp T, pari::tmp p)
{ return FpXQ_mul(y.g, x.g, T.g, p.g); }
inline pari FpXQ_pow(pari::tmp x, pari::tmp n, pari::tmp T, pari::tmp p)
{ return FpXQ_pow(x.g, n.g, T.g, p.g); }
inline pari FpXQ_powers(pari::tmp x, long l, pari::tmp T, pari::tmp p)
{ return FpXQ_powers(x.g, l, T.g, p.g); }
inline pari FpXQ_sqr(pari::tmp y, pari::tmp T, pari::tmp p)
{ return FpXQ_sqr(y.g, T.g, p.g); }
inline pari FpXQX_gcd(pari::tmp P, pari::tmp Q, pari::tmp T, pari::tmp p)
{ return FpXQX_gcd(P.g, Q.g, T.g, p.g); }
inline pari FpXQX_mul(pari::tmp x, pari::tmp y, pari::tmp T, pari::tmp p)
{ return FpXQX_mul(x.g, y.g, T.g, p.g); }
inline pari FpXQX_red(pari::tmp z, pari::tmp T, pari::tmp p)
{ return FpXQX_red(z.g, T.g, p.g); }
inline pari FpXQX_sqr(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FpXQX_sqr(x.g, T.g, p.g); }
inline pari FpXQXV_prod(pari::tmp V, pari::tmp Tp, pari::tmp p)
{ return FpXQXV_prod(V.g, Tp.g, p.g); }
inline pari FpXQYQ_pow(pari::tmp x, pari::tmp n, pari::tmp S, pari::tmp T, pari::tmp p)
{ return FpXQYQ_pow(x.g, n.g, S.g, T.g, p.g); }
inline pari FpXV_FpC_mul(pari::tmp V, pari::tmp W, pari::tmp p)
{ return FpXV_FpC_mul(V.g, W.g, p.g); }
inline pari FpXV_prod(pari::tmp V, pari::tmp p)
{ return FpXV_prod(V.g, p.g); }
inline pari FpXV_red(pari::tmp z, pari::tmp p)
{ return FpXV_red(z.g, p.g); }
inline pari FpXX_add(pari::tmp x, pari::tmp y, pari::tmp p)
{ return FpXX_add(x.g, y.g, p.g); }
inline pari FpXX_red(pari::tmp z, pari::tmp p)
{ return FpXX_red(z.g, p.g); }
inline pari FpX_rescale(pari::tmp P, pari::tmp h, pari::tmp p)
{ return FpX_rescale(P.g, h.g, p.g); }
inline pari FpY_FpXY_resultant(pari::tmp a, pari::tmp b0, pari::tmp p)
{ return FpY_FpXY_resultant(a.g, b0.g, p.g); }
inline pari Fq_inv(pari::tmp x, pari::tmp T, pari::tmp p)
{ return Fq_inv(x.g, T.g, p.g); }
inline pari Fq_invsafe(pari::tmp x, pari::tmp T, pari::tmp p)
{ return Fq_invsafe(x.g, T.g, p.g); }
inline pari Fq_mul(pari::tmp x, pari::tmp y, pari::tmp T, pari::tmp p)
{ return Fq_mul(x.g, y.g, T.g, p.g); }
inline pari Fq_neg(pari::tmp x, pari::tmp T, pari::tmp p)
{ return Fq_neg(x.g, T.g, p.g); }
inline pari Fq_neg_inv(pari::tmp x, pari::tmp T, pari::tmp p)
{ return Fq_neg_inv(x.g, T.g, p.g); }
inline pari Fq_pow(pari::tmp x, pari::tmp n, pari::tmp T, pari::tmp p)
{ return Fq_pow(x.g, n.g, T.g, p.g); }
inline pari Fq_red(pari::tmp x, pari::tmp T, pari::tmp p)
{ return Fq_red(x.g, T.g, p.g); }
inline pari FqC_to_FlxC(pari::tmp v, pari::tmp T, pari::tmp pp)
{ return FqC_to_FlxC(v.g, T.g, pp.g); }
inline pari FqM_to_FlxM(pari::tmp x, pari::tmp T, pari::tmp pp)
{ return FqM_to_FlxM(x.g, T.g, pp.g); }
inline pari FqV_roots_to_pol(pari::tmp V, pari::tmp T, pari::tmp p, long v)
{ return FqV_roots_to_pol(V.g, T.g, p.g, v); }
inline pari FqV_red(pari::tmp z, pari::tmp T, pari::tmp p)
{ return FqV_red(z.g, T.g, p.g); }
inline pari FqV_to_FlxV(pari::tmp v, pari::tmp T, pari::tmp pp)
{ return FqV_to_FlxV(v.g, T.g, pp.g); }
inline pari FqX_Fq_mul(pari::tmp P, pari::tmp U, pari::tmp T, pari::tmp p)
{ return FqX_Fq_mul(P.g, U.g, T.g, p.g); }
inline pari FqX_div(pari::tmp x, pari::tmp y, pari::tmp T, pari::tmp p)
{ return FqX_div(x.g, y.g, T.g, p.g); }
inline pari FqX_eval(pari::tmp x, pari::tmp y, pari::tmp T, pari::tmp p)
{ return FqX_eval(x.g, y.g, T.g, p.g); }
inline pari FqX_normalize(pari::tmp z, pari::tmp T, pari::tmp p)
{ return FqX_normalize(z.g, T.g, p.g); }
inline pari FqX_red(pari::tmp z, pari::tmp T, pari::tmp p)
{ return FqX_red(z.g, T.g, p.g); }
inline pari FqX_rem(pari::tmp x, pari::tmp y, pari::tmp T, pari::tmp p)
{ return FqX_rem(x.g, y.g, T.g, p.g); }
inline pari FqX_mul(pari::tmp x, pari::tmp y, pari::tmp T, pari::tmp p)
{ return FqX_mul(x.g, y.g, T.g, p.g); }
inline pari FqX_sqr(pari::tmp x, pari::tmp T, pari::tmp p)
{ return FqX_sqr(x.g, T.g, p.g); }
inline pari QXQ_inv(pari::tmp A, pari::tmp B)
{ return QXQ_inv(A.g, B.g); }
inline ulong Rg_to_Fl(pari::tmp x, ulong p)
{ return Rg_to_Fl(x.g, p); }
inline pari Rg_to_Fp(pari::tmp x, pari::tmp p)
{ return Rg_to_Fp(x.g, p.g); }
inline pari RgC_to_FpC(pari::tmp x, pari::tmp p)
{ return RgC_to_FpC(x.g, p.g); }
inline pari RgV_to_FpV(pari::tmp x, pari::tmp p)
{ return RgV_to_FpV(x.g, p.g); }
inline pari RgX_to_FpX(pari::tmp x, pari::tmp p)
{ return RgX_to_FpX(x.g, p.g); }
inline pari RgX_to_FqX(pari::tmp x, pari::tmp T, pari::tmp p)
{ return RgX_to_FqX(x.g, T.g, p.g); }
inline pari ZX_QX_resultant(pari::tmp A, pari::tmp B)
{ return ZX_QX_resultant(A.g, B.g); }
inline pari ZX_Z_add(pari::tmp y, pari::tmp x)
{ return ZX_Z_add(y.g, x.g); }
inline pari ZX_Z_mul(pari::tmp y, pari::tmp x)
{ return ZX_Z_mul(y.g, x.g); }
inline pari ZX_add(pari::tmp x, pari::tmp y)
{ return ZX_add(x.g, y.g); }
inline pari ZX_caract(pari::tmp A, pari::tmp B, long v)
{ return ZX_caract(A.g, B.g, v); }
inline pari ZX_disc(pari::tmp x)
{ return ZX_disc(x.g); }
inline int ZX_is_squarefree(pari::tmp x)
{ return ZX_is_squarefree(x.g); }
inline pari ZX_neg(pari::tmp x)
{ return ZX_neg(x.g); }
inline pari ZX_renormalize(pari::tmp x, long lx)
{ return ZX_renormalize(x.g, lx); }
inline pari ZX_resultant(pari::tmp A, pari::tmp B)
{ return ZX_resultant(A.g, B.g); }
inline pari ZX_sub(pari::tmp x, pari::tmp y)
{ return ZX_sub(x.g, y.g); }
inline pari ffinit(pari::tmp p, long n, long v)
{ return ffinit(p.g, n, v); }
inline pari from_Kronecker(pari::tmp z, pari::tmp pol)
{ return from_Kronecker(z.g, pol.g); }
inline pari init_Fq(pari::tmp p, long n, long v)
{ return init_Fq(p.g, n, v); }
inline pari modulargcd(pari::tmp a, pari::tmp b)
{ return modulargcd(a.g, b.g); }
inline pari p_stopoly(ulong m, ulong p, long v)
{ return stopoly(m, p, v); }
inline pari stopoly_gen(pari::tmp m, pari::tmp p, long v)
{ return stopoly_gen(m.g, p.g, v); }
inline pari to_Kronecker(pari::tmp P, pari::tmp Q)
{ return to_Kronecker(P.g, Q.g); }
inline int is_rational(pari::tmp x)
{ return is_rational(x.g); }
inline pari RgM_to_RgXV(pari::tmp x, long v)
{ return RgM_to_RgXV(x.g, v); }
inline pari RgM_to_RgXX(pari::tmp x, long v, long w)
{ return RgM_to_RgXX(x.g, v, w); }
inline pari RgM_zc_mul(pari::tmp x, pari::tmp y)
{ return RgM_zc_mul(x.g, y.g); }
inline pari RgM_zm_mul(pari::tmp x, pari::tmp y)
{ return RgM_zm_mul(x.g, y.g); }
inline pari RgV_to_RgX(pari::tmp x, long v)
{ return RgV_to_RgX(x.g, v); }
inline pari RgV_zc_mul(pari::tmp x, pari::tmp y)
{ return RgV_zc_mul(x.g, y.g); }
inline pari RgV_zm_mul(pari::tmp x, pari::tmp y)
{ return RgV_zm_mul(x.g, y.g); }
inline int RgX_is_rational(pari::tmp x)
{ return RgX_is_rational(x.g); }
inline pari RgX_mul(pari::tmp x, pari::tmp y)
{ return RgX_mul(x.g, y.g); }
inline pari RgX_mulspec(pari::tmp a, pari::tmp b, long na, long nb)
{ return RgX_mulspec(a.g, b.g, na, nb); }
inline pari RgX_powers(pari::tmp a, pari::tmp T, long l)
{ return RgX_powers(a.g, T.g, l); }
inline pari RgX_renormalize(pari::tmp x)
{ return RgX_renormalize(x.g); }
inline pari RgX_rescale(pari::tmp P, pari::tmp h)
{ return RgX_rescale(P.g, h.g); }
inline pari RgX_unscale(pari::tmp P, pari::tmp h)
{ return RgX_unscale(P.g, h.g); }
inline pari RgXQ_mul(pari::tmp x, pari::tmp y, pari::tmp T)
{ return RgXQ_mul(x.g, y.g, T.g); }
inline pari RgXQ_powers(pari::tmp x, long l, pari::tmp T)
{ return RgXQ_powers(x.g, l, T.g); }
inline pari RgXQ_sqr(pari::tmp x, pari::tmp T)
{ return RgXQ_sqr(x.g, T.g); }
inline pari RgXQC_red(pari::tmp P, pari::tmp T)
{ return RgXQC_red(P.g, T.g); }
inline pari RgXQV_red(pari::tmp P, pari::tmp T)
{ return RgXQV_red(P.g, T.g); }
inline pari RgXQX_RgXQ_mul(pari::tmp x, pari::tmp y, pari::tmp T)
{ return RgXQX_RgXQ_mul(x.g, y.g, T.g); }
inline pari RgXQX_mul(pari::tmp x, pari::tmp y, pari::tmp T)
{ return RgXQX_mul(x.g, y.g, T.g); }
inline pari RgXQX_red(pari::tmp P, pari::tmp T)
{ return RgXQX_red(P.g, T.g); }
inline pari RgXQX_sqr(pari::tmp x, pari::tmp T)
{ return RgXQX_sqr(x.g, T.g); }
inline pari RgXV_unscale(pari::tmp v, pari::tmp h)
{ return RgXV_unscale(v.g, h.g); }
inline pari RgX_Rg_div(pari::tmp y, pari::tmp x)
{ return RgX_Rg_div(y.g, x.g); }
inline pari RgX_Rg_mul(pari::tmp y, pari::tmp x)
{ return RgX_Rg_mul(y.g, x.g); }
inline pari RgX_RgXQ_compo(pari::tmp f, pari::tmp x, pari::tmp T)
{ return RgX_RgXQ_compo(f.g, x.g, T.g); }
inline pari RgX_mulXn(pari::tmp x, long d)
{ return RgX_mulXn(x.g, d); }
inline pari RgX_shift_shallow(pari::tmp x, long n)
{ return RgX_shift_shallow(x.g, n); }
inline pari RgX_shift(pari::tmp a, long n)
{ return RgX_shift(a.g, n); }
inline pari RgX_sqr(pari::tmp x)
{ return RgX_sqr(x.g); }
inline pari RgX_sqrspec(pari::tmp a, long na)
{ return RgX_sqrspec(a.g, na); }
inline pari RgX_to_RgV(pari::tmp x, long N)
{ return RgX_to_RgV(x.g, N); }
inline pari RgXV_to_RgM(pari::tmp v, long n)
{ return RgXV_to_RgM(v.g, n); }
inline pari RgXX_to_RgM(pari::tmp v, long n)
{ return RgXX_to_RgM(v.g, n); }
inline pari RgXY_swap(pari::tmp x, long n, long w)
{ return RgXY_swap(x.g, n, w); }
inline pari cleanroots(pari::tmp x, long l)
{ return cleanroots(x.g, l); }
inline int isrealappr(pari::tmp x, long l)
{ return isrealappr(x.g, l); }
inline pari roots(pari::tmp x, long l)
{ return roots(x.g, l); }
inline pari roots0(pari::tmp x, long flag, long l)
{ return roots0(x.g, flag, l); }
inline pari galoissubcyclo(pari::tmp N, pari::tmp sg, long flag, long v)
{ return galoissubcyclo(N.g, sg.g, flag, v); }
inline pari p_polsubcyclo(long n, long d, long v)
{ return polsubcyclo(n, d, v); }
inline pari p_subcyclo(long n, long d, long v)
{ return subcyclo(n, d, v); }
inline pari znstar_small(pari::tmp zn)
{ return znstar_small(zn.g); }
inline pari subfields(pari::tmp nf, pari::tmp d)
{ return subfields(nf.g, d.g); }
inline pari subfields0(pari::tmp nf, pari::tmp d)
{ return subfields0(nf.g, d.g); }
inline void forsubgroup(entree * oep, pari::tmp cyc, pari::tmp bound, char * och)
{ forsubgroup(oep, cyc.g, bound.g, och); }
inline pari subgrouplist(pari::tmp cyc, pari::tmp bound)
{ return subgrouplist(cyc.g, bound.g); }
inline pari bnrL1(pari::tmp bnr, pari::tmp sbgrp, long flag, long prc = prec)
{ return bnrL1(bnr.g, sbgrp.g, flag, prc); }
inline pari bnrrootnumber(pari::tmp bnr, pari::tmp chi, long flag, long prc = prec)
{ return bnrrootnumber(bnr.g, chi.g, flag, prc); }
inline pari bnrstark(pari::tmp bnr, pari::tmp subgroup, long prc = prec)
{ return bnrstark(bnr.g, subgroup.g, prc); }
inline pari direuler0(entree * ep, pari::tmp a, pari::tmp b, char * ch, pari::tmp c)
{ return direuler0(ep, a.g, b.g, ch, c.g); }
inline pari divsum(pari::tmp num, entree * ep, char * ch)
{ return divsum(num.g, ep, ch); }
inline void fordiv(pari::tmp a, entree * ep, char * ch)
{ fordiv(a.g, ep, ch); }
inline void forpari(entree * ep, pari::tmp a, pari::tmp b, char * ch)
{ forpari(ep, a.g, b.g, ch); }
inline void forprime(entree * ep, pari::tmp a, pari::tmp b, char * ch)
{ forprime(ep, a.g, b.g, ch); }
inline void forstep(entree * ep, pari::tmp a, pari::tmp b, pari::tmp s, char * ch)
{ forstep(ep, a.g, b.g, s.g, ch); }
inline void forvec(entree * ep, pari::tmp x, char * ch, long flag)
{ forvec(ep, x.g, ch, flag); }
inline pari matrice(pari::tmp nlig, pari::tmp ncol, entree * ep1, entree * ep2, char * ch)
{ return matrice(nlig.g, ncol.g, ep1, ep2, ch); }
inline pari p_polzag(long n, long m)
{ return polzag(n, m); }
inline pari p_polzagreel(long n, long m, long prc = prec)
{ return polzagreel(n, m, prc); }
inline pari prodeuler0(entree * ep, pari::tmp a, pari::tmp b, char * ch, long prc = prec)
{ return prodeuler0(ep, a.g, b.g, ch, prc); }
inline pari prodinf0(entree * ep, pari::tmp a, char * ch, long flag, long prc = prec)
{ return prodinf0(ep, a.g, ch, flag, prc); }
inline pari produit(entree * ep, pari::tmp a, pari::tmp b, char * ch, pari::tmp x)
{ return produit(ep, a.g, b.g, ch, x.g); }
inline pari somme(entree * ep, pari::tmp a, pari::tmp b, char * ch, pari::tmp x)
{ return somme(ep, a.g, b.g, ch, x.g); }
inline pari sumalt0(entree * ep, pari::tmp a, char * ch, long flag, long prc = prec)
{ return sumalt0(ep, a.g, ch, flag, prc); }
inline pari sumpos0(entree * ep, pari::tmp a, char * ch, long flag, long prc = prec)
{ return sumpos0(ep, a.g, ch, flag, prc); }
inline pari suminf0(entree * ep, pari::tmp a, char * ch, long prc = prec)
{ return suminf0(ep, a.g, ch, prc); }
inline pari vecteur(pari::tmp nmax, entree * ep, char * ch)
{ return vecteur(nmax.g, ep, ch); }
inline pari vecteursmall(pari::tmp nmax, entree * ep, char * ch)
{ return vecteursmall(nmax.g, ep, ch); }
inline pari vvecteur(pari::tmp nmax, entree * ep, char * ch)
{ return vvecteur(nmax.g, ep, ch); }
inline pari zbrent0(entree * ep, pari::tmp a, pari::tmp b, char * ch, long prc = prec)
{ return zbrent0(ep, a.g, b.g, ch, prc); }
inline pari bnfisintnorm(pari::tmp x, pari::tmp y)
{ return bnfisintnorm(x.g, y.g); }
inline pari thue(pari::tmp thueres, pari::tmp rhs, pari::tmp ne)
{ return thue(thueres.g, rhs.g, ne.g); }
inline pari thueinit(pari::tmp pol, long flag, long prc = prec)
{ return thueinit(pol.g, flag, prc); }
inline pari p_Pi2n(long n, long prc = prec)
{ return Pi2n(n, prc); }
inline pari p_PiI2(long prc = prec)
{ return PiI2(prc); }
inline pari p_PiI2n(long n, long prc = prec)
{ return PiI2n(n, prc); }
inline pari agm(pari::tmp x, pari::tmp y, long prc = prec)
{ return agm(x.g, y.g, prc); }
inline pari exp_Ir(pari::tmp x)
{ return exp_Ir(x.g); }
inline pari gcos(pari::tmp x, long prc = prec)
{ return gcos(x.g, prc); }
inline pari gcotan(pari::tmp x, long prc = prec)
{ return gcotan(x.g, prc); }
inline pari gexp(pari::tmp x, long prc = prec)
{ return gexp(x.g, prc); }
inline pari glog(pari::tmp x, long prc = prec)
{ return glog(x.g, prc); }
inline pari gpow(pari::tmp x, pari::tmp n, long prc = prec)
{ return gpow(x.g, n.g, prc); }
inline pari gpowgs(pari::tmp x, long n)
{ return gpowgs(x.g, n); }
inline pari gsin(pari::tmp x, long prc = prec)
{ return gsin(x.g, prc); }
inline pari gsqrt(pari::tmp x, long prc = prec)
{ return gsqrt(x.g, prc); }
inline pari gtan(pari::tmp x, long prc = prec)
{ return gtan(x.g, prc); }
inline pari log0(pari::tmp x, long flag, long prc = prec)
{ return log0(x.g, flag, prc); }
inline pari mpcos(pari::tmp x)
{ return mpcos(x.g); }
inline pari p_mpeuler(long prc = prec)
{ return mpeuler(prc); }
inline pari mpexp(pari::tmp x)
{ return mpexp(x.g); }
inline pari mpexp1(pari::tmp x)
{ return mpexp1(x.g); }
inline pari mplog(pari::tmp x)
{ return mplog(x.g); }
inline pari p_mplog2(long prc = prec)
{ return mplog2(prc); }
inline pari p_mppi(long prc = prec)
{ return mppi(prc); }
inline pari mpsin(pari::tmp x)
{ return mpsin(x.g); }
inline pari powiu(pari::tmp p, ulong k)
{ return powiu(p.g, k); }
inline pari p_powuu(ulong p, ulong k)
{ return powuu(p, k); }
inline pari sqrtr(pari::tmp x)
{ return sqrtr(x.g); }
inline pari sqrtnr(pari::tmp x, long n)
{ return sqrtnr(x.g, n); }
inline pari palog(pari::tmp x)
{ return palog(x.g); }
inline pari powgi(pari::tmp x, pari::tmp n)
{ return powgi(x.g, n.g); }
inline pari teich(pari::tmp x)
{ return teich(x.g); }
inline pari p_bernfrac(long n)
{ return bernfrac(n); }
inline pari p_bernreal(long n, long prc = prec)
{ return bernreal(n, prc); }
inline pari p_bernvec(long nomb)
{ return bernvec(nomb); }
inline pari gach(pari::tmp x, long prc = prec)
{ return gach(x.g, prc); }
inline pari gacos(pari::tmp x, long prc = prec)
{ return gacos(x.g, prc); }
inline pari garg(pari::tmp x, long prc = prec)
{ return garg(x.g, prc); }
inline pari gash(pari::tmp x, long prc = prec)
{ return gash(x.g, prc); }
inline pari gasin(pari::tmp x, long prc = prec)
{ return gasin(x.g, prc); }
inline pari gatan(pari::tmp x, long prc = prec)
{ return gatan(x.g, prc); }
inline pari gath(pari::tmp x, long prc = prec)
{ return gath(x.g, prc); }
inline pari gch(pari::tmp x, long prc = prec)
{ return gch(x.g, prc); }
inline pari ggamd(pari::tmp x, long prc = prec)
{ return ggamd(x.g, prc); }
inline pari ggamma(pari::tmp x, long prc = prec)
{ return ggamma(x.g, prc); }
inline pari glngamma(pari::tmp x, long prc = prec)
{ return glngamma(x.g, prc); }
inline pari gpsi(pari::tmp x, long prc = prec)
{ return gpsi(x.g, prc); }
inline pari gsh(pari::tmp x, long prc = prec)
{ return gsh(x.g, prc); }
inline pari gth(pari::tmp x, long prc = prec)
{ return gth(x.g, prc); }
inline pari p_mpfactr(long n, long prc = prec)
{ return mpfactr(n, prc); }
inline pari dilog(pari::tmp x, long prc = prec)
{ return dilog(x.g, prc); }
inline pari eint1(pari::tmp x, long prc = prec)
{ return eint1(x.g, prc); }
inline pari eta(pari::tmp x, long prc = prec)
{ return eta(x.g, prc); }
inline pari eta0(pari::tmp x, long flag, long prc = prec)
{ return eta0(x.g, flag, prc); }
inline pari gerfc(pari::tmp x, long prc = prec)
{ return gerfc(x.g, prc); }
inline pari gpolylog(long m, pari::tmp x, long prc = prec)
{ return gpolylog(m, x.g, prc); }
inline void gpolylogz(long m, pari::tmp x, pari::tmp y)
{ gpolylogz(m, x.g, y.g); }
inline pari gzeta(pari::tmp x, long prc = prec)
{ return gzeta(x.g, prc); }
inline pari hyperu(pari::tmp a, pari::tmp b, pari::tmp gx, long prc = prec)
{ return hyperu(a.g, b.g, gx.g, prc); }
inline pari incgam(pari::tmp a, pari::tmp x, long prc = prec)
{ return incgam(a.g, x.g, prc); }
inline pari incgam0(pari::tmp a, pari::tmp x, pari::tmp z, long prc = prec)
{ return incgam0(a.g, x.g, z.g, prc); }
inline pari incgam2(pari::tmp a, pari::tmp x, long prc = prec)
{ return incgam2(a.g, x.g, prc); }
inline pari incgamc(pari::tmp a, pari::tmp x, long prc = prec)
{ return incgamc(a.g, x.g, prc); }
inline pari hbessel1(pari::tmp n, pari::tmp z, long prc = prec)
{ return hbessel1(n.g, z.g, prc); }
inline pari hbessel2(pari::tmp n, pari::tmp z, long prc = prec)
{ return hbessel2(n.g, z.g, prc); }
inline pari ibessel(pari::tmp n, pari::tmp z, long prc = prec)
{ return ibessel(n.g, z.g, prc); }
inline pari jbessel(pari::tmp n, pari::tmp z, long prc = prec)
{ return jbessel(n.g, z.g, prc); }
inline pari jbesselh(pari::tmp n, pari::tmp z, long prc = prec)
{ return jbesselh(n.g, z.g, prc); }
inline pari nbessel(pari::tmp n, pari::tmp z, long prc = prec)
{ return nbessel(n.g, z.g, prc); }
inline pari jell(pari::tmp x, long prc = prec)
{ return jell(x.g, prc); }
inline pari kbessel(pari::tmp nu, pari::tmp gx, long prc = prec)
{ return kbessel(nu.g, gx.g, prc); }
inline pari kbessel0(pari::tmp nu, pari::tmp gx, long flag, long prc = prec)
{ return kbessel0(nu.g, gx.g, flag, prc); }
inline pari kbessel2(pari::tmp nu, pari::tmp x, long prc = prec)
{ return kbessel2(nu.g, x.g, prc); }
inline pari polylog(long m, pari::tmp x, long prc = prec)
{ return polylog(m, x.g, prc); }
inline pari polylog0(long m, pari::tmp x, long flag, long prc = prec)
{ return polylog0(m, x.g, flag, prc); }
inline pari polylogd(long m, pari::tmp x, long prc = prec)
{ return polylogd(m, x.g, prc); }
inline pari polylogdold(long m, pari::tmp x, long prc = prec)
{ return polylogdold(m, x.g, prc); }
inline pari polylogp(long m, pari::tmp x, long prc = prec)
{ return polylogp(m, x.g, prc); }
inline pari p_szeta(long x, long prc = prec)
{ return szeta(x, prc); }
inline pari theta(pari::tmp q, pari::tmp z, long prc = prec)
{ return theta(q.g, z.g, prc); }
inline pari thetanullk(pari::tmp q, long k, long prc = prec)
{ return thetanullk(q.g, k, prc); }
inline pari trueeta(pari::tmp x, long prc = prec)
{ return trueeta(x.g, prc); }
inline pari veceint1(pari::tmp nmax, pari::tmp C, long prc = prec)
{ return veceint1(nmax.g, C.g, prc); }
inline pari vecthetanullk(pari::tmp q, long k, long prc = prec)
{ return vecthetanullk(q.g, k, prc); }
inline pari weber0(pari::tmp x, long flag, long prc = prec)
{ return weber0(x.g, flag, prc); }
inline pari weberf(pari::tmp x, long prc = prec)
{ return weberf(x.g, prc); }
inline pari weberf1(pari::tmp x, long prc = prec)
{ return weberf1(x.g, prc); }
inline pari weberf2(pari::tmp x, long prc = prec)
{ return weberf2(x.g, prc); }
