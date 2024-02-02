#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <pthread.h>


using namespace std;

typedef unsigned char RGB[3];
using Vec3f =  parser::Vec3f ;

#define MESH 0
#define SPHERE 1
#define TRIANGLE 2

struct HitObject{
    bool hit = false;
    float t;
    int material_id;
    Vec3f normal;
    Vec3f intersected_point;
    Vec3f color; 
};

struct RayObject{
    Vec3f origin;
    Vec3f direction;
    float t;
    int depth = 0;
};

struct Args{
    parser::Camera &curr_cam;
    parser::Scene &scene;
    int &nx;
    int &ny;
    int c;
};


//For back face culling, the variable shadow_or_not is used.
//If shadow_or_not == 1, then the function is called from check_shadow function, and back face culling is not applied.

Vec3f compute_color(RayObject &r ,HitObject &hit_obj, parser::Scene &scene);


Vec3f normalization_op(Vec3f &a){
    float len = sqrt(a.x * a.x + a.y * a.y + a.z * a.z) ;
    if(len == 0) return a;
    Vec3f result{a.x / len, a.y / len, a.z / len};
    return result;
}

double find_det2(Vec3f &p1, Vec3f &p2, Vec3f &p3){
    double a,b,c,d,e,f,g,h,i;
    a = p1.x; b = p1.y; c = p1.z;
    d = p2.x; e = p2.y; f = p2.z;
    g = p3.x; h = p3.y; i = p3.z;
    double det = a * (e * i - f * h) - 
            b * (d * i - f * g) + 
            c * (d * h - e * g) ;

    return det;
}

Vec3f cross_product(Vec3f &a, Vec3f &b){
    Vec3f crossed;
    crossed.x = a.y * b.z - a.z * b.y;
    crossed.y = a.z * b.x - a.x * b.z;
    crossed.z = a.x * b.y - a.y * b.x;
    return crossed;
}

double dot_product(Vec3f &a, Vec3f &b){
    double result = a.x * b.x + a.y * b.y + a.z * b.z;
    return result;
}

Vec3f sub_vectors(Vec3f &a, Vec3f &b){
    Vec3f to_be_returned{a.x - b.x, a.y - b.y, a.z - b.z};
    
    return to_be_returned;
}

Vec3f multiply_vectors(Vec3f &a, Vec3f &b){
    Vec3f to_be_returned{a.x * b.x, a.y * b.y, a.z * b.z};
    
    return to_be_returned;
}

Vec3f scalar_multiply(Vec3f &a, float scalar){
    Vec3f to_be_returned{a.x * scalar, a.y * scalar, a.z * scalar};
    
    return to_be_returned;
}

Vec3f add_vectors(Vec3f &a, Vec3f &b){
    Vec3f to_be_returned{a.x + b.x, a.y + b.y, a.z + b.z};
    
    return to_be_returned;
}

bool isIn(float val){
    if(val <= 1 && val >= 0) return true;
    else return false;
}

// if no weight is given, give weights as 1, 1, 1
Vec3f weighted_add_3_vectors(Vec3f &vec1, Vec3f &vec2, Vec3f &vec3, float w1, float w2, float w3){
    Vec3f vec_result;
    vec_result.x = w1 * vec1.x + w2 * vec2.x + w3 * vec3.x;
    vec_result.y = w1 * vec1.y + w2 * vec2.y + w3 * vec3.y;
    vec_result.z = w1 * vec1.z + w2 * vec2.z + w3 * vec3.z;

    return vec_result;
}


bool intersection_face(parser::Face &curr_face, vector<Vec3f> &vertex_coors, Vec3f &e_cam, Vec3f &normalized_d, float &t0, float &t1,int shadow_or_not){
    
    float d_dot_n = dot_product(curr_face.normal,normalized_d);
    if(d_dot_n > 0 && !shadow_or_not) return false; 

    parser::Vec3f point_a = vertex_coors[curr_face.v0_id - 1];
    parser::Vec3f point_b = vertex_coors[curr_face.v1_id - 1];
    parser::Vec3f point_c = vertex_coors[curr_face.v2_id - 1];

    Vec3f a_minus_e = sub_vectors(point_a,e_cam);
    Vec3f a_minus_b = sub_vectors(point_a,point_b);
    Vec3f a_minus_c = sub_vectors(point_a,point_c);


    double det_A = find_det2(a_minus_b, a_minus_c, normalized_d); // |A|
    if(det_A == 0) return false; // no intersection
    double t = find_det2(a_minus_b, a_minus_c, a_minus_e) / det_A;
    if((t < t0) || (t > t1)) return false; // no intersection
    double theta = find_det2(a_minus_b, a_minus_e, normalized_d) / det_A;
    if((theta < 0) || (theta > 1)) return false; // no intersection
    double beta = find_det2(a_minus_e, a_minus_c, normalized_d) / det_A;
    if((beta < 0) || ((beta + theta) > 1)) return false; // no intersection
    
    t1 = t;
    return true;
    //double alpha = 1 - (beta + theta);
}



bool find_intersection(HitObject &hit_obj, RayObject &r,parser::Scene &scene,int shadow_or_not){
    bool hit = false;
    float t0 = 0.0;
    float t1 = INFINITY;

    int size_sphere = scene.spheres.size();
    int size_mesh = scene.meshes.size();
    int size_tri = scene.triangles.size();

    Vec3f e_cam = r.origin;
    Vec3f normalized_d = r.direction;
    std::vector<parser::Vec3f> &vertex_coors = scene.vertex_data;

    for(int m = 0; m < size_mesh; m++){
        parser::Mesh& curr_mesh = scene.meshes[m];
        int n_faces = curr_mesh.faces.size();
        
        for(int f = 0; f < n_faces; f++){
            //for each face(triangle), calculate the normal
            parser::Face& curr_face = curr_mesh.faces[f]; //current triangle

            if(intersection_face(curr_face,vertex_coors,e_cam,normalized_d,t0,t1,shadow_or_not)){
                hit_obj.hit = true;
                hit_obj.material_id = curr_mesh.material_id;
                hit_obj.t = t1;
                //hit_obj.intersected_point = weighted_add_3_vectors(point_a,point_b,point_c,alpha, beta, theta); // or e + td can be used
                Vec3f t_times_d = scalar_multiply(normalized_d,t1);
                hit_obj.intersected_point = add_vectors(e_cam,t_times_d);
                hit_obj.normal = curr_face.normal;
                hit = true;
            }                
        }

    }

    
    for(int m = 0; m < size_sphere; m++){
        parser::Sphere curr_sphere = scene.spheres[m];
        float radius = curr_sphere.radius;
        Vec3f center = vertex_coors[curr_sphere.center_vertex_id - 1];

        Vec3f e_minus_c = sub_vectors(e_cam, center);
        float d_dot_ec = dot_product(normalized_d, e_minus_c);
        float d_dot_d = dot_product(normalized_d,normalized_d);
        float ec_dot_ec = dot_product(e_minus_c,e_minus_c);



        // if delta == 0 -> corner
        // delta < 0 -> completely misses
        // d.d < 0 -> complex numbers
        float delta = sqrt(d_dot_ec * d_dot_ec - d_dot_d * (ec_dot_ec - radius * radius));
        if(d_dot_d > 0 && delta >= 0){ // intersection
            float ta = (-d_dot_ec + delta ) / d_dot_d;
            float tb = (-d_dot_ec - delta ) / d_dot_d;

            if((tb < 0) || (tb > t1)) continue; // no intersection

            
            hit_obj.hit = true;
            hit_obj.material_id = curr_sphere.material_id;
            hit_obj.t = tb;
            //hit_obj.intersected_point = Vec3f{e_cam.x + tb * normalized_d.x, e_cam.y + tb * normalized_d.y, e_cam.z + tb * normalized_d.z};
            Vec3f t_times_d = scalar_multiply(normalized_d,tb);
            hit_obj.intersected_point = add_vectors(e_cam,t_times_d);

            Vec3f p_minus_c = sub_vectors(hit_obj.intersected_point,center);
            Vec3f normal_vec = normalization_op(p_minus_c); // (p - c) / r
            hit_obj.normal = normal_vec;
            hit = true;
            t1 = tb;

            //Vec3f intersected_p1 = {e_cam.x + ta * normalized_d.x, e_cam.y + ta * normalized_d.y, e_cam.z + ta * normalized_d.z};
            //Vec3f intersected_p2 = {e_cam.x + tb * normalized_d.x, e_cam.y + tb * normalized_d.y, e_cam.z + tb * normalized_d.z};
        } 
    }



    for(int m = 0; m < size_tri; m++){
        parser::Triangle &curr_tri = scene.triangles[m];

        parser::Face &curr_face = curr_tri.indices; //current triangle face

        if(intersection_face(curr_face,vertex_coors,e_cam,normalized_d,t0,t1,0)){
            hit_obj.hit = true;
            hit_obj.material_id = curr_tri.material_id;
            hit_obj.t = t1;
            //hit_obj.intersected_point = weighted_add_3_vectors(point_a,point_b,point_c,alpha, beta, theta); // or e + td can be used
            Vec3f t_times_d = scalar_multiply(normalized_d,t1);
            hit_obj.intersected_point = add_vectors(e_cam,t_times_d);
            hit_obj.normal = curr_face.normal;
            hit = true;
        }  

    }

    r.t = t1;

    return hit;
                
}


parser::Vec3i clamp(Vec3f float_color, int min, int max){
    parser::Vec3i int_color;
    int_color.x = round(float_color.x);
    int_color.y = round(float_color.y);
    int_color.z = round(float_color.z);

    if(int_color.x < min) int_color.x = min;
    if(int_color.y < min) int_color.y = min;
    if(int_color.z < min) int_color.z = min;

    if(int_color.x > max) int_color.x = max;
    if(int_color.y > max) int_color.y = max;
    if(int_color.z > max) int_color.z = max;

    return int_color;

}

Vec3f reflect(Vec3f &w0_dir, Vec3f &normal){
    Vec3f wr;
    float n_dot_w0 = dot_product(normal,w0_dir);
    Vec3f n_with_coeff = {normal.x * 2 * n_dot_w0, normal.y * 2 * n_dot_w0, normal.z * 2 * n_dot_w0};
    wr = sub_vectors(n_with_coeff,w0_dir);
    return wr;
}

bool check_shadow(parser::Scene &scene, HitObject &hit_obj,parser::PointLight curr_light){

    //todo: send rays from intersection point to each light source
    //todo: 
    RayObject shadow_ray;
    Vec3f p = hit_obj.intersected_point;
    float epsilon = scene.shadow_ray_epsilon;
    Vec3f above_p = scalar_multiply(hit_obj.normal,epsilon);

    Vec3f light_pos = curr_light.position;
    Vec3f light_dir = sub_vectors(light_pos,p); // w_i
    //light_dir = normalization_op(light_dir); // l = w_i / |w_i|
    shadow_ray.direction = light_dir;
    shadow_ray.origin = add_vectors(p,above_p);
    HitObject shadow_hit_obj;
    if(find_intersection(shadow_hit_obj,shadow_ray,scene,1)){
        if((shadow_hit_obj.t < 1) && (shadow_hit_obj.t > 0)) return true;
    }
    return false;
}


Vec3f apply_shading(RayObject &r,HitObject &hit_obj,parser::Scene &scene){
    // intersection
    
    Vec3f e_cam = r.origin;
    Vec3f normalized_d = r.direction;
    Vec3f w0 = sub_vectors(e_cam,hit_obj.intersected_point); // w_o , below used
    w0 = normalization_op(w0);

    int material_index = hit_obj.material_id - 1;
    parser::Material curr_material = scene.materials[material_index];


    // --------------------- ambient ---------------------
    Vec3f k_ambient = curr_material.ambient;
    Vec3f ambient = scene.ambient_light;
    Vec3f La_ambient = {k_ambient.x * ambient.x, k_ambient.y * ambient.y, k_ambient.z * ambient.z};

    Vec3f color = La_ambient;

    

    if(curr_material.is_mirror == true){
        // --------------------- mirror ---------------------
        RayObject wr;
        wr.direction = reflect(w0,hit_obj.normal);
        wr.direction = normalization_op(wr.direction);
        float epsilon = scene.shadow_ray_epsilon;
        Vec3f above_p = scalar_multiply(hit_obj.normal,epsilon);
        wr.origin = add_vectors(hit_obj.intersected_point,above_p);
        //wr.origin = hit_obj.intersected_point;
        wr.depth = r.depth + 1;
        HitObject mirror_hit_obj;
        Vec3f mirrored_color = compute_color(wr,mirror_hit_obj,scene);
        mirrored_color = multiply_vectors(mirrored_color,curr_material.mirror);
        color = add_vectors(color,mirrored_color);
    }
  
    // --------------------- diffuse and specular ---------------------
    // for each light source
    int size_pointlights = scene.point_lights.size();
    for(int m = 0; m < size_pointlights; m++){
        
        parser::PointLight curr_light = scene.point_lights[m];
        
        bool isInShadow = check_shadow(scene,hit_obj,curr_light);
        if(isInShadow) continue;

        Vec3f k_diffuse = curr_material.diffuse;
        Vec3f k_specular = curr_material.specular;
        float phong_exp = curr_material.phong_exponent;

        Vec3f light_pos = curr_light.position;
        Vec3f intensity = curr_light.intensity;
        Vec3f p = hit_obj.intersected_point;
        Vec3f normal = hit_obj.normal;

        Vec3f light_dir = sub_vectors(light_pos,p); // w_i
        light_dir = normalization_op(light_dir); // l = w_i / |w_i|
        Vec3f view_dir = sub_vectors(e_cam,p); // w_o , below used
        //Vec3f view_dir = w0.direction; // w_o
        view_dir = normalization_op(view_dir);
        
        //---------------------- diffuse ----------------------
        float cos_theta = dot_product(normal,light_dir); // cos_theta_i = n.w_i
        if(cos_theta < 0) cos_theta = 0; // cos_theta_i = max(0,cos_theta_i

        Vec3f q_minus_p = sub_vectors(light_pos,p);
        float mag_qp = sqrt(q_minus_p.x * q_minus_p.x + q_minus_p.y * q_minus_p.y + q_minus_p.z * q_minus_p.z); //magnitude of q - p
        Vec3f E_d = {intensity.x / (mag_qp * mag_qp), intensity.y / (mag_qp * mag_qp), intensity.z / (mag_qp * mag_qp)}; // irradiance at distance (q-p) from the light source

        Vec3f L0_diffuse = {k_diffuse.x * E_d.x * cos_theta, k_diffuse.y * E_d.y * cos_theta, k_diffuse.z * E_d.z * cos_theta};

        //---------------------- specular ----------------------
        Vec3f half_vec = {light_dir.x + view_dir.x, light_dir.y + view_dir.y, light_dir.z + view_dir.z};
        half_vec = normalization_op(half_vec); // h = (l + w_o) / |l + w_o|
        float cos_alpha = dot_product(normal,half_vec); // cos_alpha = n.h , phong_exponent
        if(cos_alpha < 0) cos_alpha = 0; // cos_alpha = max(0,cos_alpha)

        float cos_alpha_pow = pow(cos_alpha,phong_exp);
        Vec3f L0_specular = {E_d.x * k_specular.x * cos_alpha_pow, E_d.y * k_specular.y * cos_alpha_pow, E_d.z * k_specular.z * cos_alpha_pow};
        color = {color.x + L0_diffuse.x + L0_specular.x , color.y + L0_diffuse.y + L0_specular.y , color.z + L0_diffuse.z + L0_specular.z };
        
}


    return color;
}

void find_face_normal(parser::Face &curr_face,vector<Vec3f> &vertex_coors){
    if(curr_face.ncalc == true) return;
    parser::Vec3f point_a = vertex_coors[curr_face.v0_id - 1];
    parser::Vec3f point_b = vertex_coors[curr_face.v1_id - 1];
    parser::Vec3f point_c = vertex_coors[curr_face.v2_id - 1];
    Vec3f b_minus_a = sub_vectors(point_b,point_a);
    Vec3f c_minus_a = sub_vectors(point_c,point_a);
    Vec3f normal_vec = cross_product(b_minus_a,c_minus_a);
    curr_face.normal = normalization_op(normal_vec);
    curr_face.ncalc = true;
}

Vec3f compute_color(RayObject &r ,HitObject &hit_obj, parser::Scene &scene){
    if(r.depth > scene.max_recursion_depth) return Vec3f{0,0,0};
    if(find_intersection(hit_obj,r,scene,0)) {

        return apply_shading(r,hit_obj,scene);
    }
    else if(r.depth == 0){
        //return background color
        Vec3f background_color = {(float) scene.background_color.x,(float) scene.background_color.y, (float)scene.background_color.z};
        return background_color;
    }
    else return Vec3f{0,0,0};
}

//args -> parser::Camera &curr_cam, parser::Scene &scene, int &nx, int &ny, unsigned char* image
void* render_camera(void * args){
    Args* arguments = (Args*) args; 
    
    parser::Camera &curr_cam = arguments->curr_cam;
    parser::Scene &scene = arguments->scene;
    int c = arguments->c;
    int nx = arguments->nx;
    int ny = arguments->ny;

    unsigned char* image = new unsigned char [nx * ny * 3];
    parser::Vec3f &e_cam = curr_cam.position;
    float l_cam = curr_cam.near_plane.x;
    float r_cam = curr_cam.near_plane.y;
    float b_cam = curr_cam.near_plane.z;
    float t_cam = curr_cam.near_plane.w;

    parser::Vec3f gaze = curr_cam.gaze;
    parser::Vec3f dir_v = curr_cam.up; // up
    dir_v = normalization_op(dir_v);
    parser::Vec3f dir_w = {-gaze.x, -gaze.y, -gaze.z};
    dir_w = normalization_op(dir_w);
    parser::Vec3f dir_u = cross_product(dir_v,dir_w);
    dir_u = normalization_op(dir_u);


    float &coeff_z = curr_cam.near_distance ;
    float pw = (r_cam - l_cam) / (float) nx ; // pixel width
    float ph = (t_cam - b_cam) / (float) ny ; // pixel height

    int paint_pixel = 0;
    // for each (i,j) pair such that i in [0,ny - 1] and j in [0,nx - 1]
    for(int i = 0; i < ny; i++){
        float sv = ph * (i + 0.5);
        for(int j = 0; j < nx; j++){

            HitObject hit_obj;
            RayObject r;

            float su = pw * (j + 0.5);
            
            Vec3f m_vec = {e_cam.x - dir_w.x * coeff_z, e_cam.y - dir_w.y * coeff_z, e_cam.z - dir_w.z * coeff_z};
            Vec3f q_vec = {m_vec.x + l_cam * (dir_u.x) + t_cam * dir_v.x, m_vec.y + l_cam * (dir_u.y) + t_cam * dir_v.y, m_vec.z + l_cam * (dir_u.z) + t_cam * dir_v.z};
            Vec3f s_vec = {q_vec.x + su * dir_u.x - sv * dir_v.x , q_vec.y + su * dir_u.y - sv * dir_v.y , q_vec.z + su * dir_u.z - sv * dir_v.z };
            Vec3f d_vec = sub_vectors(s_vec,e_cam);
            Vec3f normalized_d = normalization_op(d_vec);
            r.direction = normalized_d;
            r.origin = e_cam;

            Vec3f color_float = compute_color(r,hit_obj,scene);
            parser::Vec3i color = clamp(color_float,0,255);

            
            image[paint_pixel++] = color.x;
            image[paint_pixel++] = color.y;
            image[paint_pixel++] = color.z;

        }
    }
    std::string filename = scene.cameras[c].image_name;
    const char* filenameChar = filename.c_str();
    write_ppm(filenameChar, image, nx, ny);
    delete[] image;
    arguments = NULL;

    return 0;
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file

    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    std::vector<parser::Camera> &cams_vec = scene.cameras;
    std::vector<parser::Vec3f> &vertex_coors = scene.vertex_data;
    int nCameras = scene.cameras.size();
    

    // ---------------------find normals of all of the triangles---------------------
    int size_tri = scene.triangles.size();
    for(int m = 0; m < size_tri; m++){
        find_face_normal(scene.triangles[m].indices,vertex_coors);
    }

    int size_mesh = scene.meshes.size();
    for(int m = 0; m < size_mesh; m++){
        parser::Mesh& curr_mesh = scene.meshes[m];
        int n_faces = curr_mesh.faces.size();
        for(int f = 0; f < n_faces; f++){
            find_face_normal(curr_mesh.faces[f],vertex_coors);
        }
    }

    // render as many images as the number of cameras

    pthread_t threads[nCameras];
    for(int c = 0; c < nCameras; c++){
        parser::Camera &curr_cam = cams_vec[c];
        int &nx = curr_cam.image_width;
        int &ny = curr_cam.image_height;
        //unsigned char* image = new unsigned char [nx * ny * 3];

        Args* argx = new Args{curr_cam,scene,nx,ny,c};
        //Args& argx = Args{curr_cam,scene,nx,ny,image}
        pthread_create(&threads[c], NULL, render_camera, (void *) argx);

        //std::string filename = "test" + to_string(c) + ".ppm";

        
        //delete[] image;

        //delete argx;
 
    }
    for (int i = 0 ; i < nCameras; i++) 
        pthread_join(threads[i], NULL);
    
}
