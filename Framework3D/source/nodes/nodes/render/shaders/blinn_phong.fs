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

for(int i = 0; i < light_count; i ++) {

float shadow_map = texture(shadow_maps, vec3(uv, lights[i].shadow_map_id)).x;
//ambient light
float ka=0.6;
vec3 ambient=ka*lights[i].color;
//diffuse
vec3 norm=normalize(texture2D(normalMapSampler,uv).xyz);
vec3 lightDir=normalize(lights[i].position-pos);

float ks=metal*0.8;
float kd=1.0-ks;
float diff=max(dot(norm,lightDir),0.0);
vec3 diffuse=lights[i].color*(diff*kd);


vec3 viewDir=normalize(camPos-pos);
vec3 reflectDir=reflect(-lightDir,norm);
vec3 h=normalize(viewDir+lightDir);

float spec=pow(max(dot(h,norm),0.0),((1-roughness)*20));
vec3 specular =ks*lights[i].color*(spec);

//Shadow 
mat4 projection=lights[i].light_projection;
mat4 view=lights[i].light_view;
vec4 clipPos =  projection*view*(vec4(pos, 1.0));
float depth_tmp=(clipPos.z / clipPos.w);
float x_tmp=(clipPos.x/clipPos.w*0.5)+0.5;
float y_tmp=(clipPos.y/clipPos.w*0.5)+0.5;
float closestDepth = texture(shadow_maps, vec3(x_tmp,y_tmp, lights[i].shadow_map_id)).x;

float bias=0.005;
//float bias=max(0.005*(1.0-dot(lightDir,normal)),0.001);
float isS=depth_tmp-bias>closestDepth?0.0:1.0;
if(x_tmp>1.0||y_tmp>1.0||x_tmp<0||y_tmp<0)
    isS=0.0;





// PCF/PCSS
float shadow=0.0;
//float texelSize = (1.0/1024);
vec2 texelSize=  1.0/textureSize(shadow_maps,0).xy;
int filterX=2;


//PCSS part
if(depth_tmp - bias > closestDepth)
{
float d_bl=0.0;
float count_b=0.0;
for(int x = -1; x <= 1; ++x)
{
    for(int y = -1; y <= 1; ++y)
    {
        float pcfDepth = texture(shadow_maps, vec3(x_tmp+x*texelSize.x,y_tmp+y*texelSize.y, lights[i].shadow_map_id)).x;
        if(depth_tmp - bias > pcfDepth)
        {
            count_b+=1;
            d_bl+=pcfDepth;
        }
    }   
}
d_bl=d_bl/count_b;
filterX=int(25.0*d_bl/depth_tmp);
filterX=min(filterX,25);
filterX=max(filterX,0);
}
//PCSS part

for(int x = -filterX; x <= filterX; ++x)
{
    for(int y = -filterX; y <= filterX; ++y)
    {
        float pcfDepth = texture(shadow_maps, vec3(x_tmp+x*texelSize.x,y_tmp+y*texelSize.y, lights[i].shadow_map_id)).x;
        shadow += depth_tmp - bias > pcfDepth ? 0.0 : 1.0;        
    }    
}
float count=2*filterX+1;
shadow = shadow/(count*count);
if(x_tmp>1.0||y_tmp>1.0||x_tmp<0||y_tmp<0)
    shadow=0.0;
//Color=vec4(vec3(shadow),1.0);
//Color+=vec4(vec3(isS),1.0);
//isS means normal shadow,shadow means PCF shadow
vec3 result=shadow*(specular+diffuse)+ambient;
Color+=vec4(result*textureColor,1.0);




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