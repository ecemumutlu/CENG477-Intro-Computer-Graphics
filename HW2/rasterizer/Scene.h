#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Matrix4.h"


class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<std::vector<double> > depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera *camera);
	Matrix4 ModelingTransformation(Mesh meshes);
	void getCamTransformation(Camera *camera, Matrix4 &MCam);
	void getOrthoPersProjection(Camera *camera, Matrix4 &MProjection, int ProjectionType);
	void getViewportTransformation(Camera *camera, Matrix4 &MViewPort);
	Matrix4 applyTranslation(Translation *translation);
	Matrix4 applyScaling(Scaling *scaling);
	Matrix4 applyRotation(Rotation *rotation);
	bool isTriangleFacingCamera(const Vec3& v1, const Vec3& v2, const Vec3& v3, const Vec3& cameraPos);
	Vec3 crossProduct(const Vec3& v1, const Vec3& v2);
	//float dotProduct(const Vec3& v1, const Vec3& v2);
	Vec3 normalize(const Vec3& v);
	Vec3 convertVec4ToVec3(const Vec4& v);
	void rasterizeTriangle(const Triangle& triangle, int screenWidth, int screenHeight);
	Color interpolateColor(const Vec3& bc, const Triangle& triangle);
	void rasterizeTriangleSolidMode(Camera *camera, const Triangle& triangle, int screenWidth, int screenHeight);
	bool clipping(Vec4 &v0, Vec4 &v1);
	void applyPerspectiveDivide(Vec4 &point);
	void getModelingTransformation(Mesh &mesh, Matrix4 &MModel);
	void drawLine(Vec4 start, Vec4 end, Camera* camera) ;
	void setPixel(int x, int y, const Color& color, Camera* camera);
	void rasterizeTriangleSolidMode(Camera *camera, Vec3& v1, Vec3& v2, Vec3& v3, int screenWidth, int screenHeight);
  	void rasterizeWireframe(const Vec4& v1, const Vec4& v2, const Vec4& v3, const Matrix4& viewportMatrix, Camera *camera);

	Color interpolateColor2(double alpha, double beta, double gamma, Vec4 &v1, Vec4 &v2, Vec4 &v3);
	void rasterizeTriangleSolidMode2(Camera *camera, Vec4& v1, Vec4& v2, Vec4& v3, int screenWidth, int screenHeight);

	void firstCase(Vec4 startVertex, Vec4 endVertex, Vec3 startColor, Vec3 endColor);
	void secondCase(Vec4 vertexA, Vec4 vertexB, Vec3 color0, Vec3 color1);
	void thirdCase(Vec4 vertexA, Vec4 vertexB, Vec3 color0, Vec3 color1);
	void fourthCase(Vec4 vertexA, Vec4 vertexB, Vec3 color0, Vec3 color1);
};

#endif
