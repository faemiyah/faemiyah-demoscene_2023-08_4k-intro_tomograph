//#g_shader_header
#version 450

layout(location=0) uniform sampler3D u_tex;
layout(location=1) uniform sampler2D u_fbo;
layout(location=2) uniform int u_time;
#if defined(USE_LD)
layout(location=3) uniform vec3 u_control[10];
#endif

out vec4 output_color;

#if defined(USE_LD)
vec2 g_position = gl_FragCoord.xy / vec2(u_control[3].y / 2, u_control[3].z / 2) - 1;
#else
#if (DISPLAY_MODE == 1080) || (DISPLAY_MODE == -1080)
vec2 g_position = gl_FragCoord.xy / vec2(960, 540) - 1;
#else
vec2 g_position = gl_FragCoord.xy / vec2(640, 360) - 1;
#endif
#endif

#if 0

/// Box SDF.
/// See: https://iquilezles.org/articles/distfunctions/
/// \param pos Position.
/// \param bpos Box position.
/// \param sz Box size in each axis.
/// \return Signed distance.
float sdf_box(vec3 pos, vec3 bpos, vec3 sz)
{
    vec3 dist = abs(pos - bpos);
    vec3 diff = dist - sz;
    return length(max(diff,0.0)) + min(max(diff.x,max(diff.y,diff.z)),0.0);
}

/// Box frame SDF.
/// See: https://iquilezles.org/articles/distfunctions/
/// \param pos Position.
/// \param bpos Box position.
/// \param sz Box size in each axis.
/// \param fwidth Frame width.
/// \return Signed distance.
float sdf_box_frame(vec3 pos, vec3 bpos, vec3 sz, float fwidth)
{
    vec3 p = abs(pos - bpos) - sz;
    vec3 q = abs(p + fwidth) - fwidth;
    return min(min(
                length(max(vec3(p.x,q.y,q.z),0.0))+min(max(p.x,max(q.y,q.z)),0.0),
                length(max(vec3(q.x,p.y,q.z),0.0))+min(max(q.x,max(p.y,q.z)),0.0)),
            length(max(vec3(q.x,q.y,p.z),0.0))+min(max(q.x,max(q.y,p.z)),0.0));
}

/// Subtraction SDF (smooth).
/// https://iquilezles.org/articles/distfunctions/
/// \param d1 Object to subtract with.
/// \param d2 Object which is to be subtracted from.
/// \param k Smoothing radius.
/// \return Signed distance.
float sdf_smoothsubtract(float d1, float d2, float k)
{
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h);
}

/// Intersection SDF (smooth).
/// https://iquilezles.org/articles/distfunctions/
/// \param d1 Object to be intersected.
/// \param d2 Object to be intersected.
/// \param k Smoothing radius.
/// \return Signed distance.
float sdf_smoothintersect(float d1, float d2, float k)
{
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k*h*(1.0-h);
}

/// Union SDF (smooth).
/// https://iquilezles.org/articles/distfunctions/
/// \param d1 Object to be intersected.
/// \param d2 Object to be intersected.
/// \param k Smoothing radius.
/// \return Signed distance.
float sdf_smoothunion(float d1, float d2, float k)
{
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h);
}

/// Diffuse reflectance function.
/// [jboksansky]: https://github.com/boksajak/brdf/blob/master/brdf.h
/// [Oren & Nayar]: https://www1.cs.columbia.edu/CAVE/publications/pdfs/Oren_SIGGRAPH94.pdf
/// \param nor Surface normal.
/// \param ll Surface to light.
/// \param vv Surface to camera.
/// \param rsqr Roughness squared.
/// \param diffuse_reflectance Diffuse reflectance [0..1].
float brdf_diffuse(vec3 nor, vec3 ll, vec3 vv, float rsqr, float diffuse_reflectance)
{
    float ndotl = max(dot(nor, ll), 0.0);
    float ndotv = max(dot(nor, vv), 0.0);

    float A = 1.0 - 0.5 * rsqr / (rsqr + 0.33);
    float B = 0.45 * rsqr / (rsqr + 0.09);

    float anglel = acos(ndotl);
    float anglev = acos(ndotv);
    float sina = sin(max(anglel, anglev));
    float tanb = tan(min(anglel, anglev));

    vec3 projected_ll = normalize(ll - ndotl * nor);
    vec3 projected_vv = normalize(vv - ndotv * nor);
    float cos_phi_diff = dot(projected_ll, projected_vv);

    return diffuse_reflectance * ndotl * (A + B * max(0.0, cos_phi_diff) * sina * tanb);
}

/// Specular reflectance function.
/// [iY0Yi]: https://www.shadertoy.com/view/wldyRj 
/// [John Hable]: http://filmicworlds.com/blog/optimizing-ggx-shaders-with-dotlh/
/// \param nor Surface normal.
/// \param ll Surface to light.
/// \param vv Surface to camera.
/// \param rsqr Roughness squared.
/// \param F0 Fresnel term.
float brdf_specular(vec3 nor, vec3 ll, vec3 vv, float rsqr, float F0)
{
    vec3 H = normalize(vv + ll);

    float dotNL = clamp(dot(nor, ll), 0.0, 1.0);
    float dotNV = clamp(dot(nor, vv), 0.0, 1.0);
    float dotNH = clamp(dot(nor, H), 0.0, 1.0);
    float dotLH = clamp(dot(ll, H), 0.0, 1.0);

    float alphaSqr = rsqr * rsqr;
    float pi = 3.14159;
    float denom = dotNH * dotNH *(alphaSqr - 1.0) + 1.0;
    float D = alphaSqr / (pi * denom * denom);

    float dotLH5 = pow(1.0 - dotLH, 5.0);
    float F = F0 + (1.0 - F0)*(dotLH5);

    float k = rsqr * 0.5;

    return dotNL * D * F / ((dotNL * (1.0 - k) + k) * (dotNV * (1.0 - k) + k));
}

/// Reflectance function.
/// \param n Surface normal.
/// \param l Surface to light.
/// \param v Surface to camera.
/// \param roughness Roughness.
float brdf(vec3 n, vec3 l, vec3 v, float roughness)
{
    return brdf_diffuse(n, l, v, roughness, 1.0) + brdf_specular(n, l, v, roughness, 0.9);
}

#endif
