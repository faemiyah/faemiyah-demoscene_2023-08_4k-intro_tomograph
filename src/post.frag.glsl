//#g_shader_fragment_post

void main()
{
    vec2 texcoord = g_position / 2 + 0.5;
    vec3 col = texture(u_fbo, texcoord).xyz;

    // Horrible box blur.
    for(int ii = -4; ii <= 4; ii += 1)
    {
        for(int jj = -4; jj <= 4; jj += 1)
        {
            vec2 offset = vec2(ii, jj);
            vec3 tmp_col = texture(u_fbo, texcoord + offset * vec2(0.0016, 0.0028)).xyz;
            if(tmp_col.x > 8.0)
            {
                col += max((4 - length(offset)), 0) / 16 * tmp_col;
            }
        }
    }

    output_color = vec4(col / (col + 1.0), 0.7);
}
