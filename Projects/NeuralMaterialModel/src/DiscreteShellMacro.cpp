#include "../include/NeuralMaterialModel.h"
#include "../include/DiscreteShellMacro.h"

template<class NativeScaleModel>
void DiscreteShellMacro<NativeScaleModel>::addShellInplaneEnergy(T& energy)
{
    VectorXT batch_data(nFaces() * 18);
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);
        
        VectorXT nn_input(18);
        for (int i = 0; i < 3; i++)
            nn_input.segment<3>(i * 3) = vertices.row(i);
        for (int i = 0; i < 3; i++)
            nn_input.segment<3>(9 + i * 3) = undeformed_vertices.row(i);

        batch_data.segment<18>(face_idx * 18) = nn_input;
    });

    VectorXT value = native_scale_model.valueBatch(batch_data, nFaces());
    energy += value.sum();
}

inline Matrix<T, 2, 3> barycentricJacobian(const Vector<T, 3> &a, const Vector<T, 3> &b, const Vector<T, 3> &c)
{
	Vector<T, 3> v0 = b - a, v1 = c - a;
	//TV v2 = p - a;
	T d00 = v0.dot(v0); // (b - a).dot(b - a) => d d00 / d p = 0
	T d01 = v0.dot(v1); // (b - a).dot(c - a) => d d01 / d p = 0
	T d11 = v1.dot(v1); // (c - a).dot(c - a) => d d11 / dp = 0
	//T d20 = v2.dot(v0); // (p - a).dot(b - a) => d d20 / dp = (b - a)^T
	//T d21 = v2.dot(v1); //  ( p - a).dot(c - a) = > d21 / dp = (c - a)^T
	T denom = d00 * d11 - d01 * d01; // => d00 * d11 constant in => drops out, d01 constant in p => derivative is 0
	//v = (d11 * d20 - d01 * d21) / denom;
	//w = (d00 * d21 - d01 * d20) / denom;
	//u = 1.0f - v - w;
	//TV dvdp = (d11 * dd20 / dp - d01 * d d21 / dp) / denom;
	Vector<T, 3> dvdp = (d11 * (b - a) - d01 * (c - a)) / denom;
	Vector<T, 3> dwdp = (d00 * (c - a) - d01 * (b - a)) / denom;
	Matrix<T, 2, 3> result;
	result.row(0) = dvdp.transpose();
	result.row(1) = dwdp.transpose();
	return result;
}

inline Matrix<T, 3, 2> compute3DCSTDeformationGradient(
	const Vector<T, 3> &x1Undef, const Vector<T, 3> &x2Undef, const Vector<T, 3> &x3Undef,
	const Vector<T, 3> &x1, const Vector<T, 3> &x2, const Vector<T, 3> &x3)
{

	Vector<T, 3> tUndef = (x2Undef - x1Undef).normalized();
	Vector<T, 3> e2Undef = (x3Undef - x1Undef);
	Vector<T, 3> qUndef = (e2Undef - tUndef * e2Undef.dot(tUndef)).normalized();

	Matrix<T, 3, 3> x;
	x << x1, x2, x3;

	//N(b) = [1 - b1 - b2, b1, b2]
	Matrix<T, 3, 2> dNdb;
	dNdb << -1.0, -1.0,
		1.0, 0.0,
		0.0, 1.0;

	Matrix<T, 2, 3> dBdX = barycentricJacobian(x1Undef, x2Undef, x3Undef);
	Matrix<T, 3, 2> dXdXStar;
	dXdXStar << tUndef, qUndef;
	Matrix<T, 3, 2> defGrad = x * dNdb * dBdX * dXdXStar; //note that this F is not very intuitive it can contain -1 for undef configuration, but its not a problem as long as only F^T*F is used
	return defGrad;
}

template<class NativeScaleModel>
void DiscreteShellMacro<NativeScaleModel>::addShellInplaneForceEntries(VectorXT& residual)
{
    VectorXT batch_data(nFaces() * 18);
    iterateFaceParallel([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);
        FaceIdx indices = faces.segment<3>(face_idx * 3);

        T k_s = E * thickness / (1.0 - nu * nu);
        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2);
         
        Vector<T, 18> nn_input;
        for (int i = 0; i < 3; i++)
            nn_input.segment<3>(i * 3) = vertices.row(i);
        for (int i = 0; i < 3; i++)
            nn_input.segment<3>(9 + i * 3) = undeformed_vertices.row(i);
        batch_data.segment<18>(face_idx * 18) = nn_input;
    });

    VectorXT dedx = native_scale_model.gradBatch(batch_data, nFaces());

    iterateFaceSerial([&](int face_idx)
    {
        FaceIdx indices = faces.segment<3>(face_idx * 3);
        addForceEntry<3>(residual, {indices[0], indices[1], indices[2]}, -dedx.segment<9>(face_idx * 9)); 
    });
}

template<class NativeScaleModel>
void DiscreteShellMacro<NativeScaleModel>::addShellInplaneHessianEntries(std::vector<Entry>& entries)
{
    VectorXT batch_data(nFaces() * 18);
    iterateFaceParallel([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);

        FaceIdx indices = faces.segment<3>(face_idx * 3);

        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2);

        Vector<T, 18> nn_input;
        for (int i = 0; i < 3; i++)
            nn_input.segment<3>(i * 3) = vertices.row(i);
        for (int i = 0; i < 3; i++)
            nn_input.segment<3>(9 + i * 3) = undeformed_vertices.row(i);
        batch_data.segment<18>(face_idx * 18) = nn_input;
    });

    VectorXT d2edx2 = native_scale_model.hessBatch(batch_data, nFaces());

    iterateFaceSerial([&](int face_idx)
    {
        FaceIdx indices = faces.segment<3>(face_idx * 3);
        VectorXT hess_vector = d2edx2.segment<81>(face_idx * 81);
        Matrix<T, 9, 9> hess;
        for (int d = 0; d < 9; d++)
            for (int dd = 0; dd < 9; dd++)
                hess(d, dd) = hess_vector(d * 9 + dd);
        addHessianEntry<3, 3>(entries, {indices[0], indices[1], indices[2]}, hess); 
    });
}

template class DiscreteShellMacro<NeuralMaterialModel>;