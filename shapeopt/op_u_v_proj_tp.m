% OP_U_V_TP: assemble the mass matrix M = [m(i,j)], m(i,j) = (mu u_j, v_i), exploiting the tensor product structure.
%
%   mat = op_u_v_tp (spu, spv, msh, [coeff]);
%   [rows, cols, values] = op_u_v_tp (spu, spv, msh, [coeff]);
%
% INPUT:
%
%  spu:   object representing the space of trial functions (see sp_vector)
%  spv:   object representing the space of test functions (see sp_vector)
%  msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%  coeff: function handle to compute the reaction coefficient (optional)
%
% OUTPUT:
%
%  mat:    assembled mass matrix
%  rows:   row indices of the nonzero entries
%  cols:   column indices of the nonzero entries
%  values: values of the nonzero entries
% 
% Copyright (C) 2011, Carlo de Falco, Rafael Vazquez
% Copyright (C) 2016, 2017 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function varargout = op_u_v_proj_tp (space1, space2, msh, coeff)

  for icomp = 1:space1.ncomp_param
    for idim = 1:msh.ndim
      size1 = size (space1.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      size2 = size (space2.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      if (size1(2) ~= size2(2) || size1(2) ~= msh.nel_dir(idim))
        error ('One of the discrete spaces is not associated to the mesh')
      end
    end
  end

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col);
    sp2_col = sp_evaluate_col (space2, msh_col);

    if (nargin == 4)
      for idim = 1:msh.rdim
        x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
      end
      coeffs = coeff (x{:});
    else
      coeffs = ones (msh_col.nqn, msh_col.nel);
    end

    A = A + op_u_v_proj (sp1_col, sp2_col, msh_col, coeffs);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
