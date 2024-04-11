#version 430 core

// Define a uniform struct for lights
struct Light {
    // The matrices are used for shadow mapping. 
    //You need to fill it according to how we are filling it when building the normal maps (node_render_shadow_mapping.cpp). 
    // Now, they are filled with identity matrix. You need to modify C++ code in node_render_deferred_lighting.cpp.
    // Position and color are filled.
    mat4 light_projection;
    mat4 light_view;
    vec3 position;
    float radius;
    vec3 color; // Just use the same diffuse and specular color.
    int shadow_map_id;
};

layout(binding = 0) buffer lightsBuffer {
Light lights[4];
};

uniform vec2 iResolution;

uniform sampler2D diffuseColorSampler;
uniform sampler2D normalMapSampler; // You should apply normal mapping in rasterize_impl.fs
uniform sampler2D metallicRoughnessSampler;
uniform sampler2DArray shadow_maps;
uniform sampler2D position;

// uniform float alpha;
uniform vec3 camPos;

uniform int light_count;

layout(location = 0) out vec4 Color;

void main() {
vec2 uv = gl_FragCoord.xy / iResolution;

vec3 pos = texture2D(position,uv).xyz;
vec3 normal = texture2D(normalMapSampler,uv).xyz;

vec4 metalnessRoughness = texture2D(metallicRoughnessSampler,uv);
float metal = metalnessRoughness.x;
float roughness = metalnessRoughness.y;

vec3 textureColor=texture2D(diffuseColorSampler,uv).xyz;
//Color=vec4(textureColor,1.0);
//color=vec3(texture2D(diffuseColorSampler,uv).xyz);
//Color=vec4(texture2D(diffuseColorSampler,uv).xyz,1);
//color=texture2D(diffuseColorSampler,uv).xyz;
//Color=vec4(abs(texture2D(normalMapSampler,uv).xyz),1.0);
//Color=vec4(color,1.0);

//Color=vec4(ambient,1.0);
//Color=vec4((texture2D(metallicRoughnessSampler,uv).xyz),1.0);
//Color=vec4((texture2D(position,uv).xyz),1.0);

//vec3 ambient=Light.color*0.2;
//Color=vec4(ambient,1.0);

//debug
mat4 projection=lights[0].light_projection;
mat4 view=lights[0].light_view;
vec4 clipPos =  projection*(vec4(pos, 1.0));
vec3 clipT=vec3(clipPos.xyz);
//float dk=clipPos[2]/clipPos[3];
Color=vec4(clipT,1.0);
//debug
for(int i = 0; i < light_count; i ++) {

float shadow_map_value = texture(shadow_maps, vec3(uv, lights[i].shadow_map_id)).x;

float isShadowed=1.0;
vec4 clipPos = lights[i].light_projection * lights[i].light_view * vec4(pos, 1.0);
float depth=(clipPos.z / clipPos.w);
if(depth>shadow_map_value)
{
    isShadowed=0.0;
}

float ka=0.6;
vec3 ambient=ka*lights[i].color;

vec3 norm=normalize(texture2D(normalMapSampler,uv).xyz);
vec3 lightDir=normalize(lights[i].position-pos);

float kd=0.8*(roughness*2+0.5);
float diff=max(dot(norm,lightDir),0.0);
vec3 diffuse=lights[i].color*(diff*kd);

float ks=0.8;
vec3 viewDir=normalize(camPos-pos);
vec3 reflectDir=reflect(-lightDir,norm);
vec3 h=normalize(viewDir+lightDir);

float spec=pow(max(dot(h,norm),0.0),(roughness*20+10));
vec3 specular =ks*lights[i].color*(spec);

vec3 result=specular+diffuse+ambient;
//Color+=vec4(result*textureColor,1.0);



//Color=vec4(depth,0,0,1.0);
//Color=lights[i].light_projection*vec4(1.0,1.0,1.0,1.0);
//
// Visualization of shadow map
//Color = vec4(shadow_map_value,0,0,1);

// HW6_TODO: first comment the line above ("Color +=..."). That's for quick Visualization.
// You should first do the Blinn Phong shading here. You can use roughness to modify alpha. 
//Or you can pass in an alpha value through the uniform above.

// After finishing Blinn Phong shading, you can do shadow mapping with the help of the provided shadow_map_value. 
//You will need to refer to the node, node_render_shadow_mapping.cpp, for the light matrices definition. 
//Then you need to fill the mat4 light_projection; mat4 light_view; with similar approach that we fill position and color.

// For shadow mapping, as is discussed in the course, 
//you should compare the value "position depth from the light's view" against the "blocking object's depth.", 
//then you can decide whether it's shadowed.

// PCSS is also applied here.
}

}