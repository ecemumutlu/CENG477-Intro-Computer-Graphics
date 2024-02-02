#include <GL/glew.h>

class Program {
public:
    GLuint gProgram;
    GLuint vs;
    GLuint fs_normal;
    GLuint fs_ground; 

    Program() : gProgram(0), vs(0), fs_normal(0), fs_ground(0) {};
    ~Program() { };
};



