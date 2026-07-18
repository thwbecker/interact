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

#ifndef INCLUDE_POLYINTERP
#define INCLUDE_POLYINTERP

#include <vector>
#include "util/include/Matrix.hpp"

namespace util {
namespace Interp {
using namespace std;

typedef Matrix<double> Matd;
typedef vector<uint> Vecui;

// Abstract class used to wrap the user's function.
class EvalFn {
public:
  virtual ~EvalFn() {}
  // Implement this method to evaluate the user function at the n points in X,
  // where X is n x nd and F is preallocated to have n elements.
  virtual void operator()(const Matd& X, Matd& F) const
    throw (Exception) = 0;
};

class PolyMesh {
public:
  PolyMesh(const Matd& xlims, const vector<Matd>& x, const vector<Matd>& w);

  virtual ~PolyMesh() {}

  // Get an nnodes x nd array of node coordinates.
  void GetNodes(Matd& X) const;
  // Get only the new ones relative to the old PolyMesh pm0.
  void GetNewNodes(const PolyMesh& pm0, Matd& X) const;
  // Get 1-based indices to interlace old and new user function values into an
  // array of node function values for use in Eval.
  void Interlace(const PolyMesh& pm0, Vecui& Fold, Vecui& Fnew) const;

  // F is the vector of function values from EvalFnAtNodes. xi is the point at
  // which to interpolate. wrk is a work array. It's resized if necessary.
  double Eval(const Matd& F, const Matd& xi, Matd& wrk) const;
  // F has ncomp columns, where ncomp is the number of components in the
  // vector-valued function F. On return, Fe contains the interpolated value
  // for each component.
  void Eval(const Matd& F, const Matd& xi, Matd& wrk, Matd& Fe) const;

  const vector<Matd>& GetX() const { return _x; }
  const vector<Matd>& GetW() const { return _w; }
  const Matd& GetXlims() const { return _xlims; }

protected:
  // For inheriting classes:
  PolyMesh(const Matd& xlims);
  void Init();

  // nd x 2 matrix of corners of box over which interpolation is performed
  Matd _xlims;
  // Node coordinates
  vector<Matd> _x;
  // Weights
  vector<Matd> _w;
  // Some size info
  uint _nnodes; // nbr nodes
  uint _maxnn;  // max nbr nodes in a dim
};

class ChebMesh : public PolyMesh {
public:
  ChebMesh(const Matd& xlims, const Vecui& ns);
};

}}

#endif
