/* hmmvp: Software to form and apply Hierarchical Matrices
 *   Version 1.3
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 *   https://pangea.stanford.edu/research/CDFM/software
 * hmmvp is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
*/

// Andrew M. Bradley ambrad@cs.stanford.edu

#include <math.h>
#include "util/include/PolyInterp.hpp"

namespace util {
namespace Interp {

PolyMesh::
PolyMesh (const Matd& xlims, const vector<Matd>& x, const vector<Matd>& w)
  : _xlims(xlims), _x(x), _w(w)
{
  Init();
}

PolyMesh::PolyMesh (const Matd& xlims)
  : _xlims(xlims), _x(vector<Matd>(xlims.Size(1))),
    _w(vector<Matd>(xlims.Size(1)))
{}

void PolyMesh::Init () {
  size_t nd = _x.size();
  _maxnn = 0;
  _nnodes = 1;
  for (size_t i = 0; i < nd; i++) {
    size_t xs = _x[i].Size();
    if (xs > _maxnn) _maxnn = xs;
    _nnodes *= xs;
  }
}

void PolyMesh::GetNodes (Matd& X) const {
  size_t nd = _x.size();
  X.Resize(_nnodes, nd);
  for (size_t di = 0; di < nd; di++) {
    size_t rpt = 1; // nbr times to repeat _x[di](k)
    for (size_t rdi = 0; rdi < di; rdi++) rpt *= _x[rdi].Size();
    size_t k = 1, ctr = 0, nn = _x[di].Size();
    double xk = _x[di](1);
    for (size_t i = 1; i <= _nnodes; i++) {
      X(i,di+1) = xk;
      ctr++;
      if (ctr == rpt) {
        ctr = 0;
        if (k == nn) k = 1;
        else k++;
        xk = _x[di](k);
      }
    }
  }
}

void PolyMesh::GetNewNodes (const PolyMesh& pm0, Matd& X) const {
}

void PolyMesh::Interlace (const PolyMesh& pm0, Vecui& Fold, Vecui& Fnew) const {
}

// We use the barycentric formula of the second form to interpolate. See,
// for example,
//   [1] Berrut and Trefethen 2004 "Barycentric Lagrange Interpolation".
double PolyMesh::Eval (const Matd& F, const Matd& xi, Matd& wrk) const {
  // Work arrays
  wrk.Resize(_nnodes + _maxnn); // Does not realloc if the size is right
  Matd wF, ww;
  ww.SetPtr(_maxnn, 1, wrk.GetPtr());
  wF.SetPtr(_nnodes, 1, wrk.GetPtr() + _maxnn);

  size_t nd = _x.size(), Fdsz = F.Size();
  wF = F;
  double den = 1.0;
  for (size_t id = 0; id < nd; id++) { // for each dim
    size_t nn = _x[id].Size(); // nbr nodes in this dim
    size_t ei = 0; // nonzero if x[id] exactly matches an x[id]
    double swd = 0.0; // sum(wd)
    // wd = ws./(X(k,id) - xs{id});
    double xid = xi(id+1);
    for (size_t wi = 1; wi <= nn; wi++) { // for each node
      double wden = xid - _x[id](wi);
      if (wden == 0.0) { // exactly on a node?
        ei = wi;
        swd = 1.0;
        break;
      } else {
        ww(wi) = _w[id](wi) / wden;
        swd += ww(wi);
      }
    }
    den *= swd;
    size_t nch = Fdsz / nn; // nbr nn-sized chunks in Fd
    if (ei == 0) {
      // Fd(ci) <- sum(wd .* [ci'th nn-sized chunk of Fd])
      size_t os = 0;
      for (size_t ci = 0; ci < nch; ci++) {
        double wdF = 0.0; // dot(wd,[ci'th chunk of Fd])
        for (size_t wi = 1; wi <= nn; wi++) wdF += ww(wi)*wF(os+wi);
        wF(ci+1) = wdF;
        os += nn;
      }
    } else {
      // Exactly on a node in at least one dimension.
      //   Fd(ci) <- [ci'th nn-sized chunk of Fd](ei)
      size_t os = 0;
      for (size_t ci = 0; ci < nch; ci++) {
        wF(ci+1) = wF(os + ei);
        os += nn;
      }
    }
    Fdsz = nch; // new size of Fd after each chunk was collapsed
  }
  return wF(1) / den;
}

void PolyMesh::Eval (const Matd& F, const Matd& xi, Matd& wrk, Matd& Fe) const {
  size_t ncomp = F.Size(2);
  // Work arrays
  wrk.Resize(_nnodes*ncomp + _maxnn); // Does not realloc if the size is right
  Matd wF, ww;
  ww.SetPtr(_maxnn, 1, wrk.GetPtr());
  wF.SetPtr(_nnodes, ncomp, wrk.GetPtr() + _maxnn);

  size_t nd = _x.size(), Fdsz = F.Size(1);
  wF = F;
  double den = 1.0;
  for (size_t id = 0; id < nd; id++) { // for each dim
    size_t nn = _x[id].Size(); // nbr nodes in this dim
    size_t ei = 0; // nonzero if x[id] exactly matches an x[id]
    double swd = 0.0; // sum(wd)
    // wd = ws./(X(k,id) - xs{id});
    double xid = xi(id+1);
    for (size_t wi = 1; wi <= nn; wi++) { // for each node
      double wden = xid - _x[id](wi);
      if (wden == 0.0) { // exactly on a node?
        ei = wi;
        swd = 1.0;
        break;
      } else {
        ww(wi) = _w[id](wi) / wden;
        swd += ww(wi);
      }
    }
    den *= swd;
    size_t nch = Fdsz / nn; // nbr nn-sized chunks in Fd
    if (ei == 0) {
      // Fd(ci) <- sum(wd .* [ci'th nn-sized chunk of Fd])
      for (size_t cmi = 1; cmi <= ncomp; cmi++) {
        size_t os = 0;
        for (size_t ci = 0; ci < nch; ci++) {
          double wdF = 0.0; // dot(wd,[ci'th chunk of Fd])
          for (size_t wi = 1; wi <= nn; wi++) wdF += ww(wi)*wF(os+wi,cmi);
          wF(ci+1,cmi) = wdF;
          os += nn;
        }
      }
    } else {
      // Exactly on a node in at least one dimension.
      //   Fd(ci) <- [ci'th nn-sized chunk of Fd](ei)
      for (size_t cmi = 1; cmi <= ncomp; cmi++) {
        size_t os = 0;
        for (size_t ci = 0; ci < nch; ci++) {
          wF(ci+1,cmi) = wF(os + ei,cmi);
          os += nn;
        }
      }
    }
    Fdsz = nch; // new size of Fd after each chunk was collapsed
  }
  // Extract solution.
  Fe.Resize(ncomp);
  for (size_t cmi = 1; cmi <= ncomp; cmi++) Fe(cmi) = wF(1,cmi) / den;
}

// Chebyshev nodes of the second kind on [xlo xhi].
void Cheb2Nodes (double xlo, double xhi, size_t n, Matd& x, Matd& w) {
  x.Resize(n);
  w.Resize(n);
  double a = (xlo + xhi) / 2.0, b = (xhi - xlo) / 2.0;
  double sone = 1.0;
  for (size_t i = 0; i < n; i++) {
    x(n - i) = a + b*cos(double(i) * M_PI / (n - 1));
    // Barycentric weights corresponding to cheb nodes of the second kind. Not
    // normalized, but it doesn't matter when using the second form of the
    // barycentric formula.
    w(n - i) = sone;
    sone *= -1.0;
  }
  w(1) /= 2.0;
  w(n) /= 2.0;
}

ChebMesh::ChebMesh (const Matd& xlims, const Vecui& ns)
  : PolyMesh(xlims)
{
  size_t nd = ns.size();
  for (size_t id = 0; id < nd; id++)
    Cheb2Nodes(xlims(id+1,1), xlims(id+1,2), ns[id], _x[id], _w[id]);
  Init();
}

}}
