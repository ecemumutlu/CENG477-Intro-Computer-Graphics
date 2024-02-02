#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
//#include <OpenGL/gl3.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include "./sceneObjectRenderers/Object.h"
#include "./sceneObjectRenderers/Bunny.h"
#include "./sceneObjectRenderers/Quad.h"
#include "./sceneObjectRenderers/Cube.h"
#include "./sceneObjectRenderers/Texture.h"
#include "./sceneObjectRenderers/Sky.h"
#include <iostream>
#include <map>
#include <ft2build.h>
#include FT_FREETYPE_H


#define BUFFER_OFFSET(i) ((char*)NULL + (i))
#define BUNNY_YELLOW (glm::vec3(255.0f / 255.0f, 250.0f / 255.0f, 71.0f/ 255.0f))
#define CUBE_RED (glm::vec3(1.0f, 0.0f, 0.0f))
#define CUBE_YELLOW (glm::vec3(1.0f, 1.0f, 0.0f))
#define QUAD_BLUE (glm::vec3(0.0f, 0.0f, 1.0f))
#define QUAD_WHITE (glm::vec3(1.0f, 1.0f, 1.0f))
#define SCORE_YELLOW (glm::vec3(1.0f, 1.0f, 0.0f))
#define SCORE_RED (glm::vec3(1.0f, 0.0f, 0.0f))

using namespace std;

GLuint gProgram;
GLuint gProgramQuad;
GLuint gProgramText;
GLuint gProgramCube;
GLuint gProgramSky;
int gWidth, gHeight;

GLuint gScore = 0;

GLint modelingMatrixLoc[2];
GLint viewingMatrixLoc[2];
GLint projectionMatrixLoc[2];
GLint eyePosLoc[2];
GLint color[2];
GLint offsetLoc;


glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrixBunny;
glm::mat4 modelingMatrixOthers;
glm::vec3 camera;
glm::vec3 up;
glm::vec3 center;
glm::vec3 eyePos(0, 0, 0);


float scoreAcceleration = 1;
int activeProgramIndex = 1;
bool Rpressed = false;

//vector<Object> gObjects;
std::vector<Object*> gObjects;

GLuint gTextVBO;

GLint gInVertexLoc, gInNormalLoc;

//todo: erase this later, used for error
void GLAPIENTRY errorOccurredGL(GLenum source,
	GLenum type,
	GLuint id,
	GLenum severity,
	GLsizei length,
	const GLchar* message,
	const void* userParam)
{
	printf("Message from OpenGL:\nSource: 0x%x\nType: 0x%x\n"
		"Id: 0x%x\nSeverity: 0x%x\n", source, type, id, severity);
	printf("%s\n", message);
	exit(-1); // shut down the program gracefully (it does cleanup stuff too)
	// without exit(), OpenGL will constantly print the error message... make sure to kill your program.
}



/// Holds all state information relevant to a character as loaded using FreeType
struct Character {
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
    GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;


bool ParseObj(const string& fileName, glm::vec3 color, int obj_type, Object* curr_obj)
{
	fstream myfile;

	// Open the input 
	myfile.open(fileName.c_str(), std::ios::in);
	/*auto& curr_obj = gObjects.back();*/
	curr_obj->setObjectType(obj_type);

	if (myfile.is_open())
	{
		string curLine;

		while (getline(myfile, curLine))
		{
			stringstream str(curLine);
			GLfloat c1, c2, c3;
			GLuint index[9];
			string tmp;

			if (curLine.length() >= 2)
			{
				if (curLine[0] == 'v')
				{
					if (curLine[1] == 't') // texture
					{
						str >> tmp; // consume "vt"
						str >> c1 >> c2;
						curr_obj->oTextures.push_back(Texture(c1, c2));
					}
					else if (curLine[1] == 'n') // normal
					{
						str >> tmp; // consume "vn"
						str >> c1 >> c2 >> c3;
						curr_obj->oNormals.push_back(Normal(c1, c2, c3));
					}
					else // vertex
					{
						str >> tmp; // consume "v"
						str >> c1 >> c2 >> c3;
						curr_obj->oVertices.push_back(Vertex(c1 , c2 , c3 ));
						curr_obj->oColors.push_back(color);

					}
				}
				else if (curLine[0] == 'f') // face
				{
					str >> tmp; // consume "f"
					char c;
					int vIndex[3], nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0];
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1];
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2];

					assert(vIndex[0] == nIndex[0] &&
						vIndex[1] == nIndex[1] &&
						vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

					curr_obj->oFaces.push_back(Face(vIndex, tIndex, nIndex));
				}
				else
				{
					//cout << "Ignoring unidentified line in obj file: " << curLine << endl;
				}


			}

			//data += curLine;
			if (!myfile.eof())
			{
				//data += "\n";
			}
		}


		myfile.close();

		gObjects.push_back(curr_obj);
	}
	else
	{
		return false;
	}

	assert(curr_obj->oVertices.size() == curr_obj->oNormals.size());

	return true;
}


bool ReadDataFromFile(
	const string& fileName, ///< [in]  Name of the shader file
	string& data)     ///< [out] The contents of the file
{
	fstream myfile;

	// Open the input 
	myfile.open(fileName.c_str(), std::ios::in);

	if (myfile.is_open())
	{
		string curLine;

		while (getline(myfile, curLine))
		{
			data += curLine;
			if (!myfile.eof())
			{
				data += "\n";
			}
		}

		myfile.close();
	}
	else
	{
		return false;
	}

	return true;
}

GLuint createVS(const char* shaderName)
{
	string shaderSource;

	string filename(shaderName);
	if (!ReadDataFromFile(filename, shaderSource))
	{
		cout << "Cannot find file name: " + filename << endl;
		exit(-1);
	}

	GLint length = shaderSource.length();
	const GLchar* shader = (const GLchar*)shaderSource.c_str();

	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vs, 1, &shader, &length);
	glCompileShader(vs);

	char output[1024] = { 0 };
	glGetShaderInfoLog(vs, 1024, &length, output);
	//printf("VS compile log: %s\n", output);

	return vs;
}

GLuint createFS(const char* shaderName)
{
	string shaderSource;

	string filename(shaderName);
	if (!ReadDataFromFile(filename, shaderSource))
	{
		cout << "Cannot find file name: " + filename << endl;
		exit(-1);
	}

	GLint length = shaderSource.length();
	const GLchar* shader = (const GLchar*)shaderSource.c_str();

	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fs, 1, &shader, &length);
	glCompileShader(fs);

	char output[1024] = { 0 };
	glGetShaderInfoLog(fs, 1024, &length, output);
	//printf("FS compile log: %s\n", output);

	return fs;
}


void initObjectUniforms(){

	gObjects[0]->program.gProgram = gProgram;
	gObjects[0]->getUniforms();

	gObjects[1]->program.gProgram = gProgramCube;
	gObjects[1]->getUniforms();

	gObjects[2]->program.gProgram = gProgramCube;
	gObjects[2]->getUniforms();

	gObjects[3]->program.gProgram = gProgramCube;
	gObjects[3]->getUniforms();

	gObjects[4]->program.gProgram = gProgramQuad;
	gObjects[4]->getUniforms();

	gObjects[5]->program.gProgram = gProgramSky;
	gObjects[5]->getUniforms();

}

void initShaders() {
	//---------------------Normal Shader ------------------------------ //
	// Create the programs
	gProgram = glCreateProgram();

	// Create the shaders for both programs
	GLuint vs = createVS("./shaders/vert2.glsl");
	GLuint fs_normal = createFS("./shaders/frag2.glsl");

	// Attach the shaders to the programs
	glAttachShader(gProgram, vs);
	glAttachShader(gProgram, fs_normal);

	// Link the programs
	glLinkProgram(gProgram);
	GLint status;
	glGetProgramiv(gProgram, GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		std::cout << "Program link failed" << std::endl;
		exit(-1);
	}

	//--------------------------Quad Shader -------------------------- //
	gProgramQuad = glCreateProgram();
	// Create the shaders for both programs

	GLuint vs2 = createVS("./shaders/vert2.glsl");
	GLuint fs_ground = createFS("./shaders/ground_frag.glsl");


	// Attach the shaders to the programs

	glAttachShader(gProgramQuad, vs);
	glAttachShader(gProgramQuad, fs_ground);
	// Link the programs

	glLinkProgram(gProgramQuad);
	glGetProgramiv(gProgramQuad, GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		std::cout << "Program link failed" << std::endl;
		exit(-1);
	}

	//---------------------Text Shader ------------------------------ //
	// Create the programs
	gProgramText = glCreateProgram();

	// Create the shaders for both programs
	GLuint vs_text = createVS("./shaders/vert_text.glsl");
	GLuint fs_text = createFS("./shaders/frag_text.glsl");

	// Attach the shaders to the programs
	glAttachShader(gProgramText, vs_text);
	glAttachShader(gProgramText, fs_text);

	// Link the programs
	glLinkProgram(gProgramText);
	glGetProgramiv(gProgramText, GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		std::cout << "Program link failed" << std::endl;
		exit(-1);
	}

	//glBindAttribLocation(gProgramText, 3, "vertex");

	//--------------------- Cube Shader ------------------------------ //
		// Create the programs
	gProgramCube = glCreateProgram();

	// Create the shaders for both programs
	GLuint vs_cube = createVS("./shaders/vert_cube.glsl");
	GLuint fs_cube = createFS("./shaders/frag2.glsl");

	// Attach the shaders to the programs
	glAttachShader(gProgramCube, vs_cube);
	glAttachShader(gProgramCube, fs_cube);

	// Link the programs
	glLinkProgram(gProgramCube);
	glGetProgramiv(gProgramCube, GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		std::cout << "Program link failed" << std::endl;
		exit(-1);
	}

	//------------------ Texture Shader ------------------------------ //
		// Create the programs
	gProgramSky = glCreateProgram();

	// Create the shaders for both programs
	GLuint vs_texture = createVS("./shaders/vert_sky.glsl");
	GLuint fs_texture = createFS("./shaders/frag_sky.glsl");

	// Attach the shaders to the programs
	glAttachShader(gProgramSky, vs_texture);
	glAttachShader(gProgramSky, fs_texture);

	// Link the programs
	glLinkProgram(gProgramSky);
	glGetProgramiv(gProgramSky, GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		std::cout << "Program link failed" << std::endl;
		exit(-1);
	}

	initObjectUniforms();

}

void initVBOs()
{
	for (auto& obj : gObjects) {
		if (obj->object_type != SKY_TYPE){

			glGenVertexArrays(1, &obj->vao);
			assert(obj->vao > 0);
			glBindVertexArray(obj->vao);


			glEnableVertexAttribArray(0); //location 0
			glEnableVertexAttribArray(1); //location 1
			glEnableVertexAttribArray(2); //location 2

			assert(glGetError() == GL_NONE);

			//---------------

			glGenBuffers(1, &obj->gVertexAttribBuffer);
			glGenBuffers(1, &obj->gIndexBuffer);

			assert(obj->gVertexAttribBuffer > 0 && obj->gIndexBuffer > 0);

			//---------------------- Load data of Bunny first ---------------------------
			glBindBuffer(GL_ARRAY_BUFFER, obj->gVertexAttribBuffer);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, obj->gIndexBuffer);


			// The first object was bunny
			// I will have two different buffers for bunny and the others

			//Object& obj = gObjects[0];

			int oVerticesSize = obj->getVertexDataSize();
			int oNormalsSize = obj->getNormalDataSize();
			int oColorsSize = obj->getColorDataSize();
			int oFacesSize = obj->getIndexDataSize();

			int gVertexDataSizeInBytes = obj->getVertexDataSize() * 3 * sizeof(GLfloat);
			int gNormalDataSizeInBytes = obj->getNormalDataSize() * 3 * sizeof(GLfloat);
			int gColorDataSizeInBytes = obj->getColorDataSize() * 3 * sizeof(GLfloat);
			int gIndexDataSizeInBytes = obj->getIndexDataSize() * 3 * sizeof(GLuint);


			GLfloat* vertexData = new GLfloat[oVerticesSize * 3];
			GLfloat* normalData = new GLfloat[oNormalsSize * 3];
			GLfloat* colorData = new GLfloat[oColorsSize * 3];
			GLuint* indexData = new GLuint[oFacesSize * 3];


			float minX = 1e6, maxX = -1e6;
			float minY = 1e6, maxY = -1e6;
			float minZ = 1e6, maxZ = -1e6;

			vector<Vertex>& gVertices = obj->oVertices;
			vector<Normal>& gNormals = obj->oNormals;
			vector<glm::vec3>& gColors = obj->oColors;
			vector<Face>& gFaces = obj->oFaces;

			for (int i = 0; i < oVerticesSize; ++i)
			{
				vertexData[3 * i] = gVertices[i].x;
				vertexData[3 * i + 1] = gVertices[i].y;
				vertexData[3 * i + 2] = gVertices[i].z;

				minX = std::min(minX, gVertices[i].x);
				maxX = std::max(maxX, gVertices[i].x);
				minY = std::min(minY, gVertices[i].y);
				maxY = std::max(maxY, gVertices[i].y);
				minZ = std::min(minZ, gVertices[i].z);
				maxZ = std::max(maxZ, gVertices[i].z);
			}

			for (int i = 0; i < oNormalsSize; ++i)
			{
				normalData[3 * i] = gNormals[i].x;
				normalData[3 * i + 1] = gNormals[i].y;
				normalData[3 * i + 2] = gNormals[i].z;
			}

			for (int i = 0; i < oFacesSize; ++i)
			{
				indexData[3 * i] = gFaces[i].vIndex[0];
				indexData[3 * i + 1] = gFaces[i].vIndex[1];
				indexData[3 * i + 2] = gFaces[i].vIndex[2];
			}

			for (int i = 0; i < oColorsSize; ++i)
			{
				colorData[3 * i] = gColors[i].r;
				colorData[3 * i + 1] = gColors[i].g;
				colorData[3 * i + 2] = gColors[i].b;
			}


			glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes + gColorDataSizeInBytes, 0, GL_STATIC_DRAW);
			glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
			glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
			glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, gColorDataSizeInBytes, colorData);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, gIndexDataSizeInBytes, indexData, GL_STATIC_DRAW);

			// done copying to GPU memory; can free now from CPU memory
			delete[] vertexData;
			delete[] normalData;
			delete[] indexData;
			delete[] colorData;

			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
			glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes + gNormalDataSizeInBytes));

			glBindBuffer(GL_ARRAY_BUFFER, 0); //unbind
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); //unbind
		}
		else{
			glGenVertexArrays(1, &obj->vao);
			assert(obj->vao > 0);
			glBindVertexArray(obj->vao);

			glEnableVertexAttribArray(0); //location 0
			glEnableVertexAttribArray(1); //location 1

			assert(glGetError() == GL_NONE);

			//---------------

			glGenBuffers(1, &obj->gVertexAttribBuffer);
			glGenBuffers(1, &obj->gIndexBuffer);

			assert(obj->gVertexAttribBuffer > 0 && obj->gIndexBuffer > 0);

			//---------------------- Load data of Bunny first ---------------------------
			glBindBuffer(GL_ARRAY_BUFFER, obj->gVertexAttribBuffer);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, obj->gIndexBuffer);

			int oVerticesSize = obj->getVertexDataSize();
			int oFacesSize = obj->getIndexDataSize();
			int oTexturesSize = obj->getTextureDataSize();

			int gVertexDataSizeInBytes = obj->getVertexDataSize() * 3 * sizeof(GLfloat);
			int gIndexDataSizeInBytes = obj->getIndexDataSize() * 3 * sizeof(GLuint);
			int gTextureDataSizeInBytes = obj->getTextureDataSize() * 2 * sizeof(GLfloat);


			GLfloat* vertexData = new GLfloat[oVerticesSize * 3];
			GLuint* indexData = new GLuint[oFacesSize * 3];
			GLfloat* textureData = new GLfloat[oTexturesSize * 2];


			float minX = 1e6, maxX = -1e6;
			float minY = 1e6, maxY = -1e6;
			float minZ = 1e6, maxZ = -1e6;

			vector<Vertex>& gVertices = obj->oVertices;
			vector<Face>& gFaces = obj->oFaces;
			vector<Texture>& gTextures = obj->oTextures;

			for (int i = 0; i < oVerticesSize; ++i)
			{
				vertexData[3 * i] = gVertices[i].x;
				vertexData[3 * i + 1] = gVertices[i].y;
				vertexData[3 * i + 2] = gVertices[i].z;

				minX = std::min(minX, gVertices[i].x);
				maxX = std::max(maxX, gVertices[i].x);
				minY = std::min(minY, gVertices[i].y);
				maxY = std::max(maxY, gVertices[i].y);
				minZ = std::min(minZ, gVertices[i].z);
				maxZ = std::max(maxZ, gVertices[i].z);
			}

			for (int i = 0; i < oFacesSize; ++i)
			{
				indexData[3 * i] = gFaces[i].vIndex[0];
				indexData[3 * i + 1] = gFaces[i].vIndex[1];
				indexData[3 * i + 2] = gFaces[i].vIndex[2];
			}
	
			textureData[0] = gTextures[1].u;
			textureData[1] = gTextures[1].v;
			textureData[2] = gTextures[0].u;
			textureData[3] = gTextures[0].v;
			textureData[4] = gTextures[3].u;
			textureData[5] = gTextures[3].v;
			textureData[6] = gTextures[2].u;
			textureData[7] = gTextures[2].v;

			glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gTextureDataSizeInBytes , 0, GL_STATIC_DRAW);
			glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
			glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gTextureDataSizeInBytes, textureData);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, gIndexDataSizeInBytes, indexData, GL_STATIC_DRAW);

			// done copying to GPU memory; can free now from CPU memory
			delete[] vertexData;
			delete[] indexData;
			delete[] textureData;

			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
			glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

			glBindBuffer(GL_ARRAY_BUFFER, 0); //unbind
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); //unbind
		}
	}

	glGenBuffers(1, &gTextVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
	
}


void initFonts(int windowWidth, int windowHeight)
{
    // Set OpenGL options
    //glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
    glUseProgram(gProgramText);
    glUniformMatrix4fv(glGetUniformLocation(gProgramText, "projection"), 1, GL_FALSE, glm::value_ptr(projection));

    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
    }

    // Load font as face
    FT_Face face;
    if (FT_New_Face(ft, "/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf", 0, &face))
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    }

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 48);

    // Disable byte-alignment restriction
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 

    // Load first 128 characters of ASCII set
    for (GLubyte c = 0; c < 128; c++)
    {
        // Load character glyph 
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // Generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                face->glyph->bitmap.width,
                face->glyph->bitmap.rows,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                face->glyph->bitmap.buffer
                );
        // Set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // Now store character for later use
        Character character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            face->glyph->advance.x
        };
        Characters.insert(std::pair<GLchar, Character>(c, character));
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    // Destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    //
    // Configure VBO for texture quads
    //
}


void renderText(const std::string& text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color, GLuint& gScore)
{
    // Activate corresponding render state	

    glUseProgram(gProgramText);
	glActiveTexture(GL_TEXTURE0);
    glUniform3f(glGetUniformLocation(gProgramText, "textColor"), color.x, color.y, color.z);
    

	const std::string score = text + std::to_string(gScore);

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = score.begin(); c != score.end(); c++) 
    {
        Character ch = Characters[*c];

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;

        // Update VBO for each character
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },            
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }           
        };

        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID );

        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        //glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    glBindTexture(GL_TEXTURE_2D, 0);
}


void init()
{
	Bunny* bunny = new Bunny();	
	ParseObj("./sceneObjects/bunny.obj", BUNNY_YELLOW, BUNNY_TYPE, bunny);

	Cube* cube_left = new Cube(glm::vec3(-8.0f, -6.0f, -120.0f), 0); // Cube(translation, cubeID)
	ParseObj("./sceneObjects/cube.obj", CUBE_RED, CUBE_TYPE, cube_left );

	Cube* cube_middle = new Cube(glm::vec3(0.0f, -6.0f, -120.0f), 1);
	ParseObj("./sceneObjects/cube.obj", CUBE_RED, CUBE_TYPE, cube_middle );

	Cube* cube_right = new Cube(glm::vec3(8.0f, -6.0f, -120.0f), 2);
	ParseObj("./sceneObjects/cube.obj", CUBE_RED, CUBE_TYPE, cube_right);
	
	Quad* quad = new Quad(70.0f);
	ParseObj("./sceneObjects/quad.obj", QUAD_BLUE, QUAD_TYPE, quad);

	Sky* sky = new Sky("./sceneObjects/sky.jpg", 1);
	ParseObj("./sceneObjects/quad.obj", QUAD_BLUE, SKY_TYPE, sky);


	
	glEnable(GL_DEPTH_TEST);
	initShaders();
	gWidth = 1000;
	gHeight = 800;
	initFonts(gWidth, gHeight);
	initVBOs();
}




void drawModel()
{	
	glm::vec3 score_color;
	for (auto& obj : gObjects) {
		obj->drawObject();
	}
	
	//gObjects[5]->drawObject(); // draw sky first

	if(Object::redHit == true) {
		score_color = SCORE_RED;
		
	}
	else{
		score_color = SCORE_YELLOW;
	}

	renderText("Score: ", 2, 750, 1, score_color, gScore);

    assert(glGetError() == GL_NO_ERROR);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0); //unbind
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); //unbind
}


bool areFloatsEqual(float a, float b, float epsilon = 1e-4) {
	return std::abs(a - b) < epsilon;
}

void display(double deltaTime)
{

	glClearColor(0, 0, 0, 1);
	glClearDepth(1.0f);
	glClearStencil(0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	if(!Object::gameOver || Rpressed == true) {
		gObjects[0]->findRandomYellowCubeID(); // every frame, generate a random yellow cube ID but only use this random cubeID if the bunny is in the same z coordinate with the cube
		for (auto& obj : gObjects) {
			obj->setMatrices(deltaTime);
		}
		Cube::alreadyChanged = false;
		if(Object::yellowHit == true){
			gScore += 1000;
			Object::yellowHit = false;
		}	
		else if(Rpressed == true){
			Rpressed = false;
			gScore = 0;
		}	
		else{
			scoreAcceleration += CUBE_QUAD_ACCELARATION * deltaTime;
			gScore +=  std::round(scoreAcceleration);
		}

		
	}

	// Draw the scene
	drawModel();
	
}

void reshape(GLFWwindow* window, int w, int h)
{
	w = w < 1 ? 1 : w;
	h = h < 1 ? 1 : h;

	gWidth = w;
	gHeight = h;

	glViewport(0, 0, w, h);

	// Use perspective projection
	float fovyRad = (float)(90.0 / 180.0) * M_PI;
	projectionMatrix = glm::perspective(fovyRad, w / (float)h, 4.0f, 2000.0f);

	// Assume default camera position and orientation (camera is at
	// (0, 0, 0) with looking at -z direction and its up vector pointing
	// at +y direction)
	// 
	//viewingMatrix = glm::mat4(1);
	camera = glm::vec3(0, 0, 0);
	center = glm::vec3(0, 0, 0) + glm::vec3(0, 0, -1000);
	up = glm::vec3(0, 1, 0);
	
	viewingMatrix = glm::lookAt(camera, center, up);
	//viewingMatrix = glm::lookAt(glm::vec3(0, 0, 5), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));


	for (auto& obj : gObjects) {
		obj->setProjectionMatrix(projectionMatrix);
		//obj->setViewingMatrix(viewingMatrix);
	}

}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_Q && action == GLFW_PRESS)
	{
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	}
	else if (key == GLFW_KEY_R && action == GLFW_PRESS)
	{
		gScore = 0;
		scoreAcceleration = 1;
		Object::redHit = false;
		
		for(auto& obj : gObjects) {
			obj->startGame();
		}
		Object::gameOver = true;
		Rpressed = true;
	}
	else if (key == GLFW_KEY_R && action == GLFW_RELEASE)
	{
		Rpressed = false;
		Object::gameOver = false;
	}
	else if (key == GLFW_KEY_W && action == GLFW_PRESS)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else if (key == GLFW_KEY_E && action == GLFW_PRESS)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else if (key == GLFW_KEY_A && action == GLFW_PRESS)
	{
		Bunny* bunny = (Bunny*) gObjects[0];
		bunny->setPressedA(true);
	}
	else if (key == GLFW_KEY_D && action == GLFW_PRESS)
	{
		Bunny* bunny = (Bunny*) gObjects[0];
		bunny->setPressedD(true);
	}
	else if (key == GLFW_KEY_A && action == GLFW_RELEASE)
	{
		Bunny* bunny = (Bunny*) gObjects[0];
		bunny->setPressedA(false);
	}
	else if (key == GLFW_KEY_D && action == GLFW_RELEASE)
	{
		Bunny* bunny = (Bunny*) gObjects[0];
		bunny->setPressedD(false);
	}
}

void mainLoop(GLFWwindow* window)
{
	double lastTime = glfwGetTime();

	while (!glfwWindowShouldClose(window))
	{
		double currentTime = glfwGetTime(); 
		double deltaTime = currentTime - lastTime; //the result will be in seconds
		lastTime = currentTime;
		// camera.z -= 5.0f * deltaTime;
		// viewingMatrix = glm::lookAt(camera, center, up);
		display(deltaTime);
		glfwSwapBuffers(window);
		glfwPollEvents();

		

	}
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
	GLFWwindow* window;
	if (!glfwInit())
	{
		exit(-1);
	}

	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this if on MacOS

	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE); //todo: erase this later, used for error

	int width = 1000, height = 800;
	window = glfwCreateWindow(width, height, "Simple Example", NULL, NULL);

	if (!window)
	{
		glfwTerminate();
		exit(-1);
	}

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	// Initialize GLEW to setup the OpenGL Function pointers
	if (GLEW_OK != glewInit())
	{
		std::cout << "Failed to initialize GLEW" << std::endl;
		return EXIT_FAILURE;
	}



	char rendererInfo[512] = { 0 };
	strcpy(rendererInfo, (const char*)glGetString(GL_RENDERER)); // Use strcpy_s on Windows, strcpy on Linux
	strcat(rendererInfo, " - "); // Use strcpy_s on Windows, strcpy on Linux
	strcat(rendererInfo, (const char*)glGetString(GL_VERSION)); // Use strcpy_s on Windows, strcpy on Linux
	glfwSetWindowTitle(window, rendererInfo);

	init();
	

	glfwSetKeyCallback(window, keyboard);
	glfwSetWindowSizeCallback(window, reshape);

	reshape(window, width, height); // need to call this once ourselves



	// glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS); //todo: erase this later, used for error
	// glDebugMessageCallback(errorOccurredGL, NULL); //todo: erase this later, used for error

	mainLoop(window); // this does not return unless the window is closed

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
