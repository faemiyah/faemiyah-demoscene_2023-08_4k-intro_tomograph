#ifndef __g_shader_header_header__
#define __g_shader_header_header__
static const char *g_shader_header = ""
#if defined(USE_LD)
"header.glsl"
#else
"#version 450\n"
"layout(location=0)uniform sampler3D u;"
"layout(location=1)uniform sampler2D x;"
"layout(location=2)uniform int h;"
"out vec4 g;"
"vec2 s=gl_FragCoord.rg/vec2(960,540)-1;"
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
