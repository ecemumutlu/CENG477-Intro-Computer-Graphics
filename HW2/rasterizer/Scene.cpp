#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream> //BEWARE

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"
#include "Matrix4.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/

struct Vec2 {
    double x, y;
};

Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open( camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}



// Viewing Transformation Helper Functions
void Scene::getCamTransformation(Camera *camera, Matrix4 &MCam){
	// MCam-> camera transformation matrix = rotate e to align with world axes * translate e to origin

	// up vector is not always perpendicular to the gaze vector -> follow the below procedure
	camera->gaze = normalizeVec3(camera->gaze);
	camera->u = crossProductVec3(camera->gaze, camera->v);
	camera->u = normalizeVec3(camera->u);

	camera->w = inverseVec3(camera->gaze);
	camera->v = crossProductVec3(camera->u, camera->gaze);
	camera->v = normalizeVec3(camera->v);

	Vec3 cam_u = camera->u;
	Vec3 cam_v = camera->v;
	Vec3 cam_w = camera->w;
	Vec3 cam_pos = camera->position;

	//Camera transformation matrix
	MCam.values[0][0] = cam_u.x;
	MCam.values[0][1] = cam_u.y;
	MCam.values[0][2] = cam_u.z;
	MCam.values[0][3] = -1 * dotProductVec3(cam_pos, cam_u);

	MCam.values[1][0] = cam_v.x;
	MCam.values[1][1] = cam_v.y;
	MCam.values[1][2] = cam_v.z;
	MCam.values[1][3] = -1 * dotProductVec3(cam_pos, cam_v);

	MCam.values[2][0] = cam_w.x;
	MCam.values[2][1] = cam_w.y;
	MCam.values[2][2] = cam_w.z;
	MCam.values[2][3] = -1 * dotProductVec3(cam_pos, cam_w);

	MCam.values[3][0] = 0;
	MCam.values[3][1] = 0;
	MCam.values[3][2] = 0;
	MCam.values[3][3] = 1;

	return ;
}

void Scene::getOrthoPersProjection(Camera *camera, Matrix4 &MProjection, int ProjectionType){
	float left = camera->left;
	float right = camera->right;
	float bottom = camera->bottom;
	float top = camera->top;
	float near = camera->near;
	float far = camera->far;


	float r_minus_l = right - left;
	float r_plus_l = right + left;
	float t_minus_b = top - bottom;
	float t_plus_b = top + bottom;
	float f_minus_n = far - near;
	float f_plus_n = far + near;

	if(ProjectionType == ORTOGRAPHIC_PROJECTION){
		MProjection.values[0][0] = 2.0 / r_minus_l;
		MProjection.values[0][1] = 0;
		MProjection.values[0][2] = 0;
		MProjection.values[0][3] = -1.0 * r_plus_l / r_minus_l;

		MProjection.values[1][0] = 0;
		MProjection.values[1][1] = 2.0 / t_minus_b;
		MProjection.values[1][2] = 0;
		MProjection.values[1][3] = -1 * t_plus_b / t_minus_b;

		MProjection.values[2][0] = 0;
		MProjection.values[2][1] = 0;
		MProjection.values[2][2] = -2.0 / f_minus_n;
		MProjection.values[2][3] = -1.0 * f_plus_n / f_minus_n;

		MProjection.values[3][0] = 0;
		MProjection.values[3][1] = 0;
		MProjection.values[3][2] = 0;
		MProjection.values[3][3] = 1;
	}

	else if(ProjectionType == PERSPECTIVE_PROJECTION){
		MProjection.values[0][0] = 2 * near / r_minus_l;
		MProjection.values[0][1] = 0;
		MProjection.values[0][2] = r_plus_l / r_minus_l;
		MProjection.values[0][3] = 0;

		MProjection.values[1][0] = 0;
		MProjection.values[1][1] = 2 * near / t_minus_b;
		MProjection.values[1][2] = t_plus_b / t_minus_b;
		MProjection.values[1][3] = 0;

		MProjection.values[2][0] = 0;
		MProjection.values[2][1] = 0;
		MProjection.values[2][2] = -1 * (f_plus_n) / f_minus_n;
		MProjection.values[2][3] = -1 * (2 * far * near) / f_minus_n;

		MProjection.values[3][0] = 0;
		MProjection.values[3][1] = 0;
		MProjection.values[3][2] = -1;
		MProjection.values[3][3] = 0;
	}
	return;
}


void Scene::getViewportTransformation(Camera *camera, Matrix4 &MViewPort){
	int horRes = camera->horRes;
	int verRes = camera->verRes;
	
	MViewPort.values[0][0] = horRes / 2.0;
	MViewPort.values[0][1] = 0;
	MViewPort.values[0][2] = 0;
	MViewPort.values[0][3] = (horRes - 1) / 2.0;

	MViewPort.values[1][0] = 0;
	MViewPort.values[1][1] = verRes / 2.0;
	MViewPort.values[1][2] = 0;
	MViewPort.values[1][3] = (verRes - 1) / 2.0;

	MViewPort.values[2][0] = 0;
	MViewPort.values[2][1] = 0;
	MViewPort.values[2][2] = 0.5;
	MViewPort.values[2][3] = 0.5;

	MViewPort.values[3][0] = 0;
	MViewPort.values[3][1] = 0;
	MViewPort.values[3][2] = 0;
	MViewPort.values[3][3] = 1;

	return ;
}

void Scene::applyPerspectiveDivide(Vec4 &point){
	if(point.t != 1 || point.t != 0){//apply perspective divide
		float w = point.t;
		point.x = point.x / w;
		point.y = point.y / w;
		point.z = point.z / w;
		point.t = 1.0;
	}
}

bool isVisible(double den, double num, double tE, double tL){
	if(den > 0){ // potentially entering
		double t = num / den;
		if(t > tL){
			return false;
		}
		else if(t > tE){
			tE = t;
		}
	}
	else if(den < 0){ // potentially leaving
		double t = num / den;
		if(t < tE){
			return false;
		}
		else if(t < tL){
			tL = t;
		}
	}
	else if(num > 0){ // parallel and outside
		return false;
	}

	return true;
}

bool Scene::clipping(Vec4 &v0, Vec4 &v1){
	// Liang-Barsky Algorithm
	// given the line v0, v1 and the clipping window (xmin, xmax, ymin, ymax, zmin, zmax)
	// returns false if the line is outside of the clipping window
	// returns true if the line is inside of the clipping window

	// screen coordinates
	double xmin = -1; double xmax = 1;
	double ymin = -1; double ymax = 1;
	double zmin = -1; double zmax = 1;

	double tE = 0; //entering
	double tL = 1; // leaving

	bool visible = false;

	double dx = v1.x - v0.x;
	double dy = v1.y - v0.y;
	double dz = v1.z - v0.z;

	if(isVisible(dx, xmin - v0.x, tE, tL)){
		if(isVisible(-dx, v0.x - xmax, tE, tL)){
			if(isVisible(dy, ymin - v0.y, tE, tL)){
				if(isVisible(-dy, v0.y - ymax, tE, tL)){
					if(isVisible(dz, zmin - v0.z, tE, tL)){
						if(isVisible(-dz, v0.z - zmax, tE, tL)){
							visible = true;
							if(tL < 1){
								v1.x = v0.x + tL * dx;
								v1.y = v0.y + tL * dy;
								v1.z = v0.z + tL * dz;
							}
							if(tE > 0){
								v0.x = v0.x + tE * dx;
								v0.y = v0.y + tE * dy;
								v0.z = v0.z + tE * dz;
							}
						}
					}
				}
			}
		}
	}

	

	return visible;
}

Matrix4 MatrixInverse(Matrix4 &matrix){
	//orthonormal = transpose

	Matrix4 to_be_returned;
	to_be_returned.values[0][0] = matrix.values[0][0];
	to_be_returned.values[0][1] = matrix.values[1][0];
	to_be_returned.values[0][2] = matrix.values[2][0];
	to_be_returned.values[0][3] = matrix.values[3][0];

	to_be_returned.values[1][0] = matrix.values[0][1];
	to_be_returned.values[1][1] = matrix.values[1][1];
	to_be_returned.values[1][2] = matrix.values[2][1];
	to_be_returned.values[1][3] = matrix.values[3][1];

	to_be_returned.values[2][0] = matrix.values[0][2];
	to_be_returned.values[2][1] = matrix.values[1][2];
	to_be_returned.values[2][2] = matrix.values[2][2];
	to_be_returned.values[2][3] = matrix.values[3][2];

	to_be_returned.values[3][0] = matrix.values[0][3];
	to_be_returned.values[3][1] = matrix.values[1][3];
	to_be_returned.values[3][2] = matrix.values[2][3];
	to_be_returned.values[3][3] = matrix.values[3][3];

	return to_be_returned;

}

void getScalingMatrix(Matrix4 &tmpMatrix, Scaling *scaling){
    tmpMatrix = getIdentityMatrix(); //Initialize to Identity Matrix
    tmpMatrix.values[0][0] = scaling->sx; // Scaling by x-axis
    tmpMatrix.values[1][1] = scaling->sy; // Scaling by y-axis
    tmpMatrix.values[2][2] = scaling->sz; // Scaling by z-axis

}

void getTransformationMatrix(Matrix4 &tmpMatrix, Translation* translation){
    tmpMatrix = getIdentityMatrix(); //Initialize to Identity Matrix
    tmpMatrix.values[0][3] = translation->tx; // Translation along x-axis
    tmpMatrix.values[1][3] = translation->ty; // Translation along y-axis
    tmpMatrix.values[2][3] = translation->tz; // Translation along z-axis
}

void getRotationMatrix(Matrix4 &tmpMatrix, Rotation* rotation){

	Vec3 u_vec(rotation->ux, rotation->uy, rotation->uz);
	u_vec = normalizeVec3(u_vec);
	// Convert angle to radians
    double rad = (rotation->angle * (M_PI)) / 180.0;
    double cosA = cos(rad);
    double sinA = sin(rad);

	// ----------------- get v and w vectors ----------------
	int smallestIndex = 0;
	double min_val = abs(u_vec.x);

	if(abs(u_vec.y) < min_val) {
		smallestIndex = 1;
		min_val = abs(u_vec.y);
	}
	if(abs(u_vec.z) < min_val) {
		smallestIndex = 2;
		min_val = abs(u_vec.z);
	}

	Vec3 v_vec(0, 0, 0);
	v_vec.changeNthComponent(smallestIndex, 0);
	v_vec.changeNthComponent((smallestIndex + 1) % 3, -1 * u_vec.getNthComponent((smallestIndex + 2) % 3));
	v_vec.changeNthComponent((smallestIndex + 2) % 3, u_vec.getNthComponent((smallestIndex + 1) % 3));


	Vec3 w_vec = crossProductVec3(u_vec, v_vec);
	u_vec = normalizeVec3(u_vec);
	v_vec = normalizeVec3(v_vec);
	w_vec = normalizeVec3(w_vec);
	

	// ------------------ rotate uvw to align with xyz --------------
	Matrix4 alignMatrix;
	alignMatrix.values[0][0] = u_vec.x;
	alignMatrix.values[0][1] = u_vec.y;
	alignMatrix.values[0][2] = u_vec.z;
	alignMatrix.values[0][3] = 0;
	alignMatrix.values[1][0] = v_vec.x;
	alignMatrix.values[1][1] = v_vec.y;
	alignMatrix.values[1][2] = v_vec.z;
	alignMatrix.values[1][3] = 0;
	alignMatrix.values[2][0] = w_vec.x;
	alignMatrix.values[2][1] = w_vec.y;
	alignMatrix.values[2][2] = w_vec.z;
	alignMatrix.values[2][3] = 0;
	alignMatrix.values[3][0] = 0;
	alignMatrix.values[3][1] = 0;
	alignMatrix.values[3][2] = 0;
	alignMatrix.values[3][3] = 1;

	// ------------------ rotate around x axis ----------------------

	Matrix4 rotationMatrix;
	rotationMatrix.values[0][0] = 1;
	rotationMatrix.values[0][1] = 0;
	rotationMatrix.values[0][2] = 0;
	rotationMatrix.values[0][3] = 0;
	rotationMatrix.values[1][0] = 0;
	rotationMatrix.values[1][1] = cosA;
	rotationMatrix.values[1][2] = -1 * sinA;
	rotationMatrix.values[1][3] = 0;
	rotationMatrix.values[2][0] = 0;
	rotationMatrix.values[2][1] = sinA;
	rotationMatrix.values[2][2] = cosA;
	rotationMatrix.values[2][3] = 0;
	rotationMatrix.values[3][0] = 0;
	rotationMatrix.values[3][1] = 0;
	rotationMatrix.values[3][2] = 0;
	rotationMatrix.values[3][3] = 1;

	// ------------------ rotate back to uvw ------------------------
	Matrix4 alignMatrixInverse = MatrixInverse(alignMatrix);

	tmpMatrix = multiplyMatrixWithMatrix(alignMatrixInverse, rotationMatrix);
	tmpMatrix = multiplyMatrixWithMatrix(tmpMatrix, alignMatrix);

}

//ecem
void Scene::getModelingTransformation(Mesh &mesh, Matrix4 &MModel){
	MModel = getIdentityMatrix();
	for(int i = 0; i < mesh.numberOfTransformations; i++){
		char transformationType = mesh.transformationTypes[i];
		int transformationId = mesh.transformationIds[i] - 1;
		Matrix4 tmpMatrix;
		if(transformationType == 's'){
			getScalingMatrix(tmpMatrix, this->scalings[transformationId]);
		}
		else if(transformationType == 't'){
			getTransformationMatrix(tmpMatrix, this->translations[transformationId]);
		}
		else if(transformationType == 'r'){
			getRotationMatrix(tmpMatrix, this->rotations[transformationId]);
		}

		MModel = multiplyMatrixWithMatrix(tmpMatrix, MModel);
	}
}




//Backface Culling Helper Function -- DEBUG
bool Scene::isTriangleFacingCamera(const Vec3& v1, const Vec3& v2, const Vec3& v3, const Vec3& cameraPos) {
    // Calculate the normal of the triangle
	Vec3 edge1 = subtractVec3(v2, v1);
    Vec3 edge2 = subtractVec3(v3, v1);
    Vec3 normal = crossProductVec3(edge1, edge2);
    normal = normalizeVec3(normal); // Ensure the normal is a unit vector

	//Vec3 center = {(v1.x + v2.x + v3.x) / 3.0, (v1.y + v2.y + v3.y) / 3.0, (v1.z + v2.z + v3.z) / 3.0};
    float dot = dotProductVec3(normal, v1);

    // If dot product is positive, the triangle is facing the camera
    return dot > 0;
}

//Rasterization Helper Function
void Scene::setPixel(int x, int y, const Color& color, Camera* camera) {
    if (x >= 0 && x < camera->horRes && y >= 0 && y < camera->verRes) {
        this->image[x][y].r = std::max(0.0, std::min(255.0, round(color.r)));
        this->image[x][y].g = std::max(0.0, std::min(255.0, round(color.g)));
        this->image[x][y].b = std::max(0.0, std::min(255.0, round(color.b)));
    }
}


//Solid Mode Rasterization Helper Functions
void findBoundingBox2(Vec3 &v1, Vec3& v2, Vec3 &v3, Vec2& bboxmin, Vec2& bboxmax) {
	bboxmin.x = min(min(v1.x, v2.x), v3.x) ;
	bboxmin.y = min(min(v1.y, v2.y), v3.y) ;
	bboxmax.x = max(max(v1.x, v2.x), v3.x);
	bboxmax.y = max(max(v1.y, v2.y), v3.y);

}


double calculateEdge(Vec4 &v0, Vec4& v1, int x, int y){
	double to_be = x * (v0.y - v1.y) + y * (v1.x - v0 .x) + v0.x * v1.y - v1.x * v0.y;
	//double to_be_returned = x * (v1.y - v2.y) + y * (v2.x - v1.x) + v1.x * v2.y - v2.x * v1.y;
	return to_be;

}

void calculateBarycentric(Vec4& v0, Vec4& v1, Vec4& v2, int x, int y, double& alpha, double& beta, double& gamma){
	// alpha = calculateEdge(v1, v2, x, y) / calculateEdge(v1, v2, v0.x, v0.y);
	// beta = calculateEdge(v2, v0, x, y) / calculateEdge(v2, v0, v1.x, v1.y);
	// gamma = calculateEdge(v0, v1, x, y) / calculateEdge(v0, v1, v2.x, v2.y);

	float detT = (v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y);
    alpha = ((v1.y - v2.y) * (x - v2.x) + (v2.x - v1.x) * (y - v2.y)) / detT;
    beta = ((v2.y - v0.y) * (x - v2.x) + (v0.x - v2.x) * (y - v2.y)) / detT;
    gamma = 1 - alpha - beta;
}

Color Scene::interpolateColor(const Vec3& bc, const Triangle& triangle) {
    // Retrieve the color for each vertex using the colorId
    const Color& color1 = *colorsOfVertices[vertices[triangle.vertexIds[0]]->colorId];
    const Color& color2 = *colorsOfVertices[vertices[triangle.vertexIds[1]]->colorId];
    const Color& color3 = *colorsOfVertices[vertices[triangle.vertexIds[2]]->colorId];

    // Interpolate the color using barycentric coordinates
    double r = bc.x * color1.r + bc.y * color2.r + bc.z * color3.r;
    double g = bc.x * color1.g + bc.y * color2.g + bc.z * color3.g;
    double b = bc.x * color1.b + bc.y * color2.b + bc.z * color3.b;

    // Clamp the values to the range [0, 255] and round them
    r = std::max(0.0, std::min(255.0, round(r)));
    g = std::max(0.0, std::min(255.0, round(g)));
    b = std::max(0.0, std::min(255.0, round(b)));

    // Create a new Color object with the interpolated values
    return Color(r, g, b);
}

Color Scene::interpolateColor2(double alpha, double beta, double gamma, Vec4 &v1, Vec4 &v2, Vec4 &v3) {


    // Retrieve the color for each vertex using the colorId
    Color color1 = *colorsOfVertices[v1.colorId - 1];
    Color color2 = *colorsOfVertices[v2.colorId - 1];
    Color color3 = *colorsOfVertices[v3.colorId - 1];
    // Interpolate the color using barycentric coordinates
    double r = alpha * color1.r + beta * color2.r + gamma * color3.r;
    double g = alpha * color1.g + beta * color2.g + gamma * color3.g;
    double b = alpha * color1.b + beta * color2.b + gamma * color3.b;

    // Clamp the values to the range [0, 255] and round them
    r = std::max(0.0, std::min(255.0, round(r)));
    g = std::max(0.0, std::min(255.0, round(g)));
    b = std::max(0.0, std::min(255.0, round(b)));

    // Create a new Color object with the interpolated values
    return Color(r, g, b);
}





void Scene::rasterizeTriangleSolidMode2(Camera *camera, Vec4& v1, Vec4& v2, Vec4& v3, int screenWidth, int screenHeight) {
	//edge line equations

	double verRes = camera->verRes;
	double horRes = camera->horRes;

	double max_v1x = max(v1.x,0.);
	double max_v2x = max(v2.x,0.);
	double max_v3x = max(v3.x,0.);
	double min_v1x = min(v1.x,horRes);
	double min_v2x = min(v2.x,horRes);
	double min_v3x = min(v3.x,horRes);

	double max_v1y = max(v1.y,0.);
	double max_v2y = max(v2.y,0.);
	double max_v3y = max(v3.y,0.);
	double min_v1y = min(v1.y,verRes);
	double min_v2y = min(v2.y,verRes);
	double min_v3y = min(v3.y,verRes);

	double xmin = min(min(max_v1x, max_v2x), max_v3x);
	double xmax = max(max(min_v1x, min_v2x), min_v3x);

	double ymin = min(min(max_v1y, max_v2y), max_v3y) ;
	double ymax = max(max(min_v1y, min_v2y), min_v3y);


	// double xmin = min(min(v1.x, v2.x), v3.x);
	// double xmax = max(max(v1.x, v2.x), v3.x);

	// double ymin = min(min(v1.y, v2.y), v3.y);
	// double ymax = max(max(v1.y, v2.y), v3.y);


	for (int y = ymin; y <= ymax && y >= 0 && y < verRes; y++){
		for (int x = xmin; x <= xmax && x >= 0 && x < horRes; x++){
			double alpha, beta, gamma;
			calculateBarycentric(v1, v2, v3, x, y, alpha, beta, gamma);
			if (alpha >= 0 && beta >= 0 && gamma >= 0) {
				// Interpolate the color using barycentric coordinates				
				Color color = interpolateColor2(alpha, beta, gamma, v1, v2, v3); // Modify to use the colors from the vertices
				double curr_depth = alpha * v1.z + beta * v2.z + gamma * v3.z;
				//cout << curr_depth << endl;
				//assignColorToPixel(x, y, color);
				if (curr_depth < this->depth[x][y] && curr_depth >= 0) {
					//cout << "inside" << endl;
				 	this->depth[x][y] = curr_depth;
					assignColorToPixel(x, y, color);
					//cout << "still inside" << endl;
				}

			}
		}
	}

	// cout << "here" << endl;
	// for (int y = 0; y < verRes; y++){
	// 	for (int x = 0; x < horRes; x++){
	// 		assignColorToPixel(x, y, this->colorsOfVertices[x][y]);		
	// 	}
	// }
						
}






//Wireframe Rasterization Helper Function
void Scene::drawLine(Vec4 start, Vec4 end, Camera* camera) {
    const Color& startColor = *this->colorsOfVertices[start.colorId - 1];
    const Color& endColor = *this->colorsOfVertices[end.colorId - 1];

    int x0 = static_cast<int>(start.x);
    int y0 = static_cast<int>(start.y);
    int x1 = static_cast<int>(end.x);
    int y1 = static_cast<int>(end.y);

    int dx = abs(x1 - x0);
	int sx = x0 < x1 ? 1 : -1; //1 if right, -1 if left
    int dy = -abs(y1 - y0);
	int sy = y0 < y1 ? 1 : -1; //1 if down, -1 if up
    int err = (dx > dy ? dx : -dy) / 2, e2;

    while (true) {
        // Calculate the interpolation factor t
        double t = sqrt(pow(x0 - start.x, 2) + pow(y0 - start.y, 2)) / sqrt(pow(x1 - start.x, 2) + pow(y1 - start.y, 2));
        double r = (1 - t) * startColor.r + t * endColor.r;
        double g = (1 - t) * startColor.g + t * endColor.g;
        double b = (1 - t) * startColor.b + t * endColor.b;

        // Set the pixel color using setPixel
        this->setPixel(x0, y0, Color(r, g, b), camera);

        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 > -dx) { err -= dy; x0 += sx; }
        if (e2 < dy) { err += dx; y0 += sy; }
    }
}


void Scene::firstCase(Vec4 vertexA, Vec4 vertexB, Vec3 color0, Vec3 color1) {
	int vertexAyInt= static_cast<int>(vertexA.y);
	int vertexByInt= static_cast<int>(vertexB.y);
	int vertexAxInt= static_cast<int>(vertexA.x);
	int vertexBxInt= static_cast<int>(vertexB.x);
	int ByMinusAy = static_cast<int>(vertexB.y - vertexA.y);
	int BxMinusAx = static_cast<int>(vertexB.x - vertexA.x);

    int minY = std::min(vertexAyInt, vertexByInt);
	int maxX = std::max(vertexAxInt, vertexBxInt);
    int theDifference = 2 * std::abs(ByMinusAy) - std::abs(BxMinusAx);

    for (int x = std::min(vertexAxInt, vertexBxInt); x <= maxX; x++) {
        image[x][minY].g = static_cast<unsigned char>(std::round(color0.y * std::abs(vertexB.x - x) / std::abs(vertexB.x - vertexA.x) + color1.y * std::abs(vertexA.x - x) / std::abs(vertexB.x - vertexA.x)));
        image[x][minY].b = static_cast<unsigned char>(std::round(color0.z * std::abs(vertexB.x - x) / std::abs(vertexB.x - vertexA.x) + color1.z * std::abs(vertexA.x - x) / std::abs(vertexB.x - vertexA.x)));
		image[x][minY].r = static_cast<unsigned char>(std::round(color0.x * std::abs(vertexB.x - x) / std::abs(vertexB.x - vertexA.x) + color1.x * std::abs(vertexA.x - x) / std::abs(vertexB.x - vertexA.x)));
        
		if (theDifference > 0) {
			minY++; 
			theDifference += 2 * std::abs(static_cast<int>(vertexB.y - vertexA.y)) - 2 * std::abs(static_cast<int>(vertexB.x - vertexA.x));
        } else {
            theDifference += 2 * std::abs(static_cast<int>(vertexB.y - vertexA.y));
        }
    }
}
 //

void Scene::secondCase(Vec4 vertexA, Vec4 vertexB, Vec3 color0, Vec3 color1) {
	int vertexAyInt= static_cast<int>(vertexA.y);
	int vertexByInt= static_cast<int>(vertexB.y);
	int vertexAxInt= static_cast<int>(vertexA.x);
	int vertexBxInt= static_cast<int>(vertexB.x);
	int ByMinusAy = static_cast<int>(vertexB.y - vertexA.y);
	int BxMinusAx = static_cast<int>(vertexB.x - vertexA.x);

	int maxY = std::max(vertexAyInt, vertexByInt);
    int minX = std::min(vertexAxInt, vertexBxInt);
    int theDifference = 2 * std::abs(BxMinusAx) - std::abs(ByMinusAy);
    
    for (int y = std::min(vertexAyInt, vertexByInt); y <= maxY; y++) {
        image[minX][y].r = static_cast<unsigned char>(std::round(color0.x * std::abs(vertexB.y - y) / std::abs(vertexB.y - vertexA.y) + color1.x * std::abs(vertexA.y - y) / std::abs(vertexB.y - vertexA.y)));
        image[minX][y].g = static_cast<unsigned char>(std::round(color0.y * std::abs(vertexB.y - y) / std::abs(vertexB.y - vertexA.y) + color1.y * std::abs(vertexA.y - y) / std::abs(vertexB.y - vertexA.y)));
        image[minX][y].b = static_cast<unsigned char>(std::round(color0.z * std::abs(vertexB.y - y) / std::abs(vertexB.y - vertexA.y) + color1.z * std::abs(vertexA.y - y) / std::abs(vertexB.y - vertexA.y)));

        if (theDifference > 0) {
			minX++;
			theDifference += 2 * std::abs(static_cast<int>(vertexB.x - vertexA.x)) - 2 * std::abs(static_cast<int>(vertexB.y - vertexA.y));
        } else {
            theDifference += 2 * std::abs(static_cast<int>(vertexB.x - vertexA.x));
        }
    }
}

void Scene::thirdCase(Vec4 vertexA, Vec4 vertexB, Vec3 color0, Vec3 color1) {
	int vertexAyInt= static_cast<int>(vertexA.y);
	int vertexByInt= static_cast<int>(vertexB.y);
	int vertexAxInt= static_cast<int>(vertexA.x);
	int vertexBxInt= static_cast<int>(vertexB.x);
	int ByMinusAy = static_cast<int>(vertexB.y - vertexA.y);
	int BxMinusAx = static_cast<int>(vertexB.x - vertexA.x);

    int maxY = std::max(vertexAyInt, vertexByInt);
	int maxX = std::max(vertexAxInt, vertexBxInt);
    int theDifference = 2 * std::abs(ByMinusAy) - std::abs(BxMinusAx);
    
    for (int x = std::min(vertexAxInt, vertexBxInt); x <= maxX; x++) {
        image[x][maxY].r = static_cast<unsigned char>(std::round(color0.x * std::abs(vertexB.x - x) / std::abs(vertexB.x - vertexA.x) + color1.x * std::abs(vertexA.x - x) / std::abs(vertexB.x - vertexA.x)));
        image[x][maxY].g = static_cast<unsigned char>(std::round(color0.y * std::abs(vertexB.x - x) / std::abs(vertexB.x - vertexA.x) + color1.y * std::abs(vertexA.x - x) / std::abs(vertexB.x - vertexA.x)));
        image[x][maxY].b = static_cast<unsigned char>(std::round(color0.z * std::abs(vertexB.x - x) / std::abs(vertexB.x - vertexA.x) + color1.z * std::abs(vertexA.x - x) / std::abs(vertexB.x - vertexA.x)));

        if (theDifference > 0) {
			maxY--;
			theDifference += 2 * std::abs(static_cast<int>(vertexB.y - vertexA.y)) - 2 * std::abs(BxMinusAx);
        } else {
            theDifference += 2 * std::abs(ByMinusAy);
        }
    }
}

void Scene::fourthCase(Vec4 vertexA, Vec4 vertexB, Vec3 color0, Vec3 color1) {
    int vertexAyInt= static_cast<int>(vertexA.y);
	int vertexByInt= static_cast<int>(vertexB.y);
	int vertexAxInt= static_cast<int>(vertexA.x);
	int vertexBxInt= static_cast<int>(vertexB.x);
	int ByMinusAy = static_cast<int>(vertexB.y - vertexA.y);
	int BxMinusAx = static_cast<int>(vertexB.x - vertexA.x);

	int maxY = std::max(vertexAyInt, vertexByInt);
	int maxX = std::max(vertexAxInt, vertexBxInt);
    int theDifference = 2 * std::abs(BxMinusAx) - std::abs(ByMinusAy);
    
    for (int y = std::min(static_cast<int>(vertexA.y), static_cast<int>(vertexB.y)); y <= maxY; y++) {
        image[maxX][y].r = static_cast<unsigned char>(std::round(color0.x * std::abs(vertexB.y - y) / std::abs(vertexB.y - vertexA.y) + color1.x * std::abs(vertexA.y - y) / std::abs(vertexB.y - vertexA.y)));
        image[maxX][y].g = static_cast<unsigned char>(std::round(color0.y * std::abs(vertexB.y - y) / std::abs(vertexB.y - vertexA.y) + color1.y * std::abs(vertexA.y - y) / std::abs(vertexB.y - vertexA.y)));
        image[maxX][y].b = static_cast<unsigned char>(std::round(color0.z * std::abs(vertexB.y - y) / std::abs(vertexB.y - vertexA.y) + color1.z * std::abs(vertexA.y - y) / std::abs(vertexB.y - vertexA.y)));

        if (theDifference > 0) {
			maxX--;
			theDifference += 2 * std::abs(static_cast<int>(vertexB.x - vertexA.x)) - 2 * std::abs(static_cast<int>(vertexB.y - vertexA.y));
        } else {
            theDifference += 2 * std::abs(static_cast<int>(vertexB.x - vertexA.x));
        }
    }
}

//Wireframe Rasterization Function
void Scene::rasterizeWireframe(const Vec4& v1, const Vec4& v2, const Vec4& v3, const Matrix4& viewportMatrix, Camera *camera) {
    Vec4 vertices[] = {v1, v2, v3};

    for (int i = 0; i < 2; i++) {
        for (int j = i + 1; j < 3; j++) {
            Vec4 vertexA = vertices[i];
            Vec4 vertexB = vertices[j];

			Vec3 color0 = Vec3(this->colorsOfVertices[vertexA.colorId - 1]->r,
                                           this->colorsOfVertices[vertexA.colorId - 1]->g,
                                           this->colorsOfVertices[vertexA.colorId - 1]->b, -1);
            Vec3 color1 = Vec3(this->colorsOfVertices[vertexB.colorId - 1]->r,
                                           this->colorsOfVertices[vertexB.colorId - 1]->g,
                                           this->colorsOfVertices[vertexB.colorId - 1]->b, -1);
                         

            // Apply clipping
            if (!clipping(vertexA, vertexB)) {
                continue; // Skip this line if it's outside the clipping window
            }
			
			// Apply viewport transformation
            vertexA = multiplyMatrixWithVec4(viewportMatrix, vertexA);
            vertexB = multiplyMatrixWithVec4(viewportMatrix, vertexB);



						double slope = (vertexB.y - vertexA.y) / (vertexB.x - vertexA.x);
                        if (slope > 0 && slope < 1){
							firstCase(vertexA, vertexB, color0, color1);
                        }
                        else if (slope > 1){
							secondCase(vertexA, vertexB, color0, color1);
                            
                        }
                        else if (slope <= 0 && slope >= -1){
							thirdCase(vertexA, vertexB, color0, color1);
                            
                        }
                        else if (slope < -1) {
							fourthCase(vertexA, vertexB, color0, color1);
							
                        }

            // Draw the line between the vertices
			//drawLine(vertexA, vertexB, camera); // Update to pass camera
        }
    }
}




//Forward Rendering Pipeline Function
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// ########################## VIEWING TRANSFORMATIONS ##########################

	//------------ camera transformation ----------------
	Matrix4 MCam;
	this->getCamTransformation(camera, MCam);

	// --------------- projection transformations -------------
	Matrix4 MProjection;
	this->getOrthoPersProjection(camera, MProjection, camera->projectionType);

	// ---------- viewport transformation -----------------
	Matrix4 MViewPort;
	this->getViewportTransformation(camera, MViewPort);


	//##############################################################################

	Matrix4 MProCam = multiplyMatrixWithMatrix(MProjection, MCam);
	int meshesNumber = static_cast<int>(this->meshes.size());
	for(int i=0; i<meshesNumber; i++) {
		Mesh &curr_mesh = *(this->meshes)[i];
		int trianglesNumber = static_cast<int>(curr_mesh.triangles.size());

		
		// --------------- modeling transformation begin-----------------
		Matrix4 MModel;
		this->getModelingTransformation(curr_mesh, MModel);

		// --------------- modeling transformation end -----------------
		Matrix4 MProCamModel = multiplyMatrixWithMatrix(MProCam, MModel);
		for(int j=0; j<trianglesNumber; j++) {
			Triangle triangle = curr_mesh.triangles[j];


			int v1_index = triangle.vertexIds[0] - 1;
			int v2_index = triangle.vertexIds[1] - 1;
			int v3_index = triangle.vertexIds[2] - 1;

			Color v1c = *(this->colorsOfVertices[v1_index]);
			Color v2c = *(this->colorsOfVertices[v2_index]);
			Color v3c = *(this->colorsOfVertices[v3_index]);

			Vec3 *firstVertexPointer = this->vertices[v1_index];
			Vec3 *secondVertexPointer = this->vertices[v2_index];
			Vec3 *thirtVertexPointer = this->vertices[v3_index];

            Vec4 firstVertex = Vec4(firstVertexPointer->x, firstVertexPointer->y, firstVertexPointer->z, 1, firstVertexPointer->colorId);
            Vec4 secondVertex = Vec4(secondVertexPointer->x, secondVertexPointer->y, secondVertexPointer->z, 1, secondVertexPointer->colorId);
            Vec4 thirdVertex = Vec4(thirtVertexPointer->x, thirtVertexPointer->y, thirtVertexPointer->z, 1, thirtVertexPointer->colorId);
            
			//MULTIPLY

			
			firstVertex = Vec4(multiplyMatrixWithVec4(MProCamModel, firstVertex));
			secondVertex = Vec4(multiplyMatrixWithVec4(MProCamModel, secondVertex));
			thirdVertex = Vec4(multiplyMatrixWithVec4(MProCamModel, thirdVertex));

			// ----------- perspective division -------------------
			if(camera->projectionType == PERSPECTIVE_PROJECTION){
				// apply MPers * MCam * MModel to points
				// Lets say our point is p = (x, y, z, 1)
				applyPerspectiveDivide(firstVertex);
				applyPerspectiveDivide(secondVertex);
				applyPerspectiveDivide(thirdVertex);
			}

			//todo: do clipping if wireframe, dont do it if solid
			//todo: clipping color should change
			// ----------- clipping -----------------------------
			// Liang-Barsky Algorithm
			if (((this->meshes)[i])->type == WIREFRAME_MESH) { //WIREFRAME
				bool clip_line1 = this->clipping(firstVertex, secondVertex);
				bool clip_line2 = this->clipping(secondVertex, thirdVertex);
				bool clip_line3 = this->clipping(thirdVertex, firstVertex);
			}
			// // Directly convert transformed Vec4 vertices to Vec3

			Vec3 v1(firstVertex.x, firstVertex.y, firstVertex.z);
        	Vec3 v2(secondVertex.x, secondVertex.y, secondVertex.z);
        	Vec3 v3(thirdVertex.x, thirdVertex.y, thirdVertex.z);


			// ----------- culling ------------------------------
            if (this->cullingEnabled == 1 && !isTriangleFacingCamera(v1, v2, v3, camera->position)) {
                continue; // Skip rendering this triangle
            }


			// ----------- rasterization -------------------------
						//rasterizeImage
			if (((this->meshes)[i])->type == WIREFRAME_MESH) { //WIREFRAME
				rasterizeWireframe(firstVertex, secondVertex, thirdVertex, MViewPort, camera);
            }	

			if (((this->meshes)[i])->type == SOLID_MESH) { //SOLID
							// ----------- viewport transformation ---------------
			// apply MViewPort to points
				firstVertex = Vec4(multiplyMatrixWithVec4(MViewPort, firstVertex));
				secondVertex = Vec4(multiplyMatrixWithVec4(MViewPort, secondVertex));
				thirdVertex = Vec4(multiplyMatrixWithVec4(MViewPort, thirdVertex));
				rasterizeTriangleSolidMode2(camera, firstVertex, secondVertex, thirdVertex, camera->horRes, camera->verRes);

			}
		}

	}

}










