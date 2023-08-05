#ifndef __g_shader_fragment_post_header__
#define __g_shader_fragment_post_header__
static const char *g_shader_fragment_post = ""
#if defined(USE_LD)
"post.frag.glsl"
#else
"void main()"
"{"
"vec2 t=s/2+.5;"
"vec3 e=texture(x,t).rgb;"
"for(int c=-4;"
"c<=4;"
"c+=1)for(int v=-4;"
"v<=4;"
"v+=1)"
"{"
"vec2 a=vec2(c,v);"
"vec3 c=texture(x,t+a*vec2(.0016,.0028)).rgb;"
"if(c.r>8.)e+=max(4-length(a),0)/16*c;"
"}"
"g=vec4(e/(e+1.),.7);"
"}"
#endif
"";
#if !defined(DNLOAD_RENAME_UNUSED)
#if defined(__GNUC__)
#define DNLOAD_RENAME_UNUSED __attribute__((unused))
#else
#define DNLOAD_RENAME_UNUSED
#endif
#endif
#endif
