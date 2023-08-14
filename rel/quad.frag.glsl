//#g_shader_fragment_quad

mat3 rot = mat3(0.94, -0.27, 0.2, -0.18, -0.9, -0.41, 0.29, 0.35, -0.89);
int g_object_mask = 3;

/// World description.
///
/// Format is triplets of vec3s.
/// 0: position
/// 1: size (distances to edges)
/// 2: settings
///
/// Setting X: type
/// 1: Inverted quad.
/// 2: Regular quad.
/// 4: Frame.
/// 8: Light.
/// 16: Scan.
vec3 g_world[] =
{
    vec3(0), vec3(2, 1, 2),
    vec3(0, -0.4, 0), vec3(0.2, 0.05, 0.8),
    vec3(0, -0.1, -0.6), vec3(0.44, 0.44, 0.11),
    vec3(-0.5, 0.8, -0.8), vec3(0.02, 0.01, 0.06),
    vec3(0.5, 0.8, 0.8), vec3(0.02, 0.01, 0.06),
    vec3(-0.5, 0.8, -0), vec3(0.02, 0.01, 0.06),
    vec3(-0.5, 0.8, 0.8), vec3(0.02, 0.01, 0.06),
    vec3(0.5, 0.8, -0.8), vec3(0.02, 0.01, 0.06),
    vec3(0.5, 0.8, 0), vec3(0.02, 0.01, 0.06),
    vec3(0, -0.1, -0.6), vec3(0.4, 0.4, 0.02),
    vec3(0, -0.1, 0), vec3(0.4, 0.4, 0.7)
};

/// Camera data.
/// 0: Duration, frame position start, frame position end.
/// 1: Position start.
/// 2: Lookat start.
/// 3: Position end.
/// 4: Lookat end.
vec3 g_camera[] =
{
    vec3(14, -0.6, -0.6), vec3(-1.2, 0.1, 1.8), vec3(0.7, 0, -0.8), vec3(-1.5, 0, 1.5), vec3(0.7, 0, -0.8),
    vec3(10, 0.6, -0.6), vec3(-1.5, 0.1, 0.8), vec3(0.9, -0.2, -0.5), vec3(-1.7, 0.1, 0.3), vec3(1, -0.2, -0.4),
    vec3(12, -0.6, -0.6), vec3(1.4, -0.7, 1.3), vec3(-0.6, 0.3, -0.9), vec3(1.3, -0.7, 1.4), vec3(-0.5, 0.3, -0.9),
    vec3(4, -0.6, -0.4), vec3(1.3, -0.7, 1.4), vec3(-0.5, 0.3, -0.9), vec3(1.3, -0.5, 1.4), vec3(-0.5, 0.2, -0.9),
    vec3(12, -0.4, 0.6), vec3(0, -0.1, 0.1), vec3(0, 0, -1), vec3(0, -0.1, 1.2), vec3(0, 0, -1),
    vec3(10, 0.6, -0.6), vec3(0.9, 0.7, -1.4), vec3(-0.8, -0.5, 0.9), vec3(1.2, 0.5, -1.4), vec3(-0.7, -0.4, 0.7),
    vec3(16, -0.6, 0.6), vec3(0.9, -0.1, 0.3), vec3(-0.6, 0, -0.9), vec3(0.8, -0.1, 1.5), vec3(-0.6, 0, -0.9),
    vec3(14, 0.6, -0.6), vec3(-1.6, 0.2, 1.4), vec3(0.8, -0.2, -0.7), vec3(-1.9, 0.2, 0.5), vec3(1, -0.1, -0.5),
    vec3(10, -0.6, -0.6), vec3(-1.5, -0.7, 0.6), vec3(0.9, 0.4, -0.5), vec3(-1.4, -0.7, 0.7), vec3(0.8, 0.4, -0.6)
};

/// Masking data.
/// 0: Duration.
/// 1: OR mask constant.
/// 2: OR mask start.
/// 3: OR mask end.
/// 4: Shift left, negative for u_time masking off.
int g_masking[] =
{
    2, 0, 0, 0, 0,
    11, 0, 0, 63, 3,
    4, 504, 0, 0, 0,
    4, 504, 0, 63, 1,
    6, 511, 0, 0, 0,
    3, 4, 63, 0, -1,
    1, 4, 0, 0, 0,
    1, 4, 0, 1, 9,
    21, 516, 0, 0, 0,
    8, 4, 0, 63, 9,
    3, 516, 0, 0, 0,
    12, 4, 0, 63, 9,
    3, 516, 0, 0, 0,
    4, 516, 0, 63, 3,
    6, 511, 63, 0, -9,
    7, 511, 0, 0, 0,
    3, 4, 63, 0, 0,
    3, 0, 0, 0, 0
};

/// Frame size for the frame thing.
float i_frame_width = 0.033;

/// Number of concrete elements.
int i_world_elements_concrete = 3;

/// Number of total elements.
int i_world_elements_total = 11;

/// Scan index.
int i_world_element_index_tardigrade = i_world_elements_total - 1;

/// Scan index.
int i_world_element_index_scan = i_world_elements_total - 2;

/// Pi constant.
float g_pi = 3.14159;

/// Ellipsoid SDF.
/// https://iquilezles.org/articles/ellipsoids/
/// \param pos Position relative to ellipsoid.
/// \param r Ellipsoid radii.
/// \return Signed distance.
float sdf_ellipsoid(vec3 pos, vec3 r)
{
    float kk = length(pos / r);
    return kk * (kk - 1) / length(pos / (r * r));
}

/// Tardigrade SDF.
/// Distance estimator from http://www.imajeenyus.com/mathematics/20121112_distance_estimates/distance_estimation_method_for_fractals.pdf
/// \param pos Position relative to tardigrade.
float sdf_tardigrade(vec3 pos)
{
    int i_MAXITER = 8;
    float i_BAILOUT = 2;
    float i_XYSCALE = 5;
    int i_ZSCALE = 12;
    vec2 wormslices[] =     // Julia slice coordinates
    {
        vec2(0.1, 0.2),     // paeae loppu
        vec2(0.05, 0.15),   // paeae
        vec2(0.1, 0.4),     // hartia 1
        vec2(0.15, 0),      // jalka 1
        vec2(0.05, 0.45),   // hartia 2
        vec2(0.2, 0),       // jalka 2
        vec2(0, 0.5),       // hartia 3, PUOLIVAELI
        vec2(0.2, 0),       // jalka 3
        vec2(0.05, 0.45),   // hartia 4
        vec2(0.15, 0),      // jalka 4
        vec2(0.1, 0.4),     // hartia 5
        vec2(0.15, 0),      // perse
        vec2(0.3, 0),       // perse loppu
    };

    float zpos = (pos.z + 0.5) * i_ZSCALE;
    vec2 wormslice;
    float enddistance = 0;

    if(zpos < 0)
    {
        wormslice = wormslices[0];
        enddistance = -zpos;
    }
    else if(zpos > i_ZSCALE)
    {
        wormslice = wormslices[i_ZSCALE];
        enddistance = zpos - i_ZSCALE;
    }
    else {
        float mixconst = sin((zpos - floor(zpos)) * g_pi / 2);
        wormslice = mix(wormslices[int(floor(zpos))], wormslices[int(ceil(zpos))], mixconst);
    }

    float wormsegmentscale = i_XYSCALE / (pow(sin(pos.z * 1.3 - 1.6), 3) + 0.05);

    vec3 tpos = pos + vec3(0, 0.27 * sin(pos.z * 2 - g_pi / 2) + 0.24, 0.5); // Translate to center, then translate as a function of position
    //vec3 tpos = (pos + vec3(0.0, 0.27*sin(pos.z*2-1.57)+0.22, 0.5) - bpos); // Translate to center, then translate as a function of position

    vec2 z1 = vec2(tpos.x, tpos.y) * wormsegmentscale;
    vec2 z2 = vec2(tpos.x, tpos.y) * wormsegmentscale;
    vec2 dz1 = vec2(1, 0);
    vec2 dz2 = dz1;

    for(int ii=0; ii<i_MAXITER; ii++) {
        vec2 z1_new = vec2(z1.x * z1.x - z1.y * z1.y, z1.x * z1.y + z1.y * z1.x) + vec2(wormslice.x, wormslice.y);
        vec2 dz1_new = vec2(z1.x * dz1.x - z1.y * dz1.y, z1.x * dz1.y + z1.y * dz1.x) * 2;
        vec2 z2_new = vec2(z2.x * z2.x - z2.y * z2.y, z2.x * z2.y + z2.y * z2.x) + vec2(wormslice.x, -wormslice.y);
        vec2 dz2_new = vec2(z2.x * dz2.x - z2.y * dz2.y, z2.x * dz2.y + z2.y * dz2.x) * 2;
        z1 = z1_new;
        dz1 = dz1_new;
        z2 = z2_new;
        dz2 = dz2_new;

        // Symmetric with respect to Y -> same bailout for both
        if(max(length(z1), length(z2)) > i_BAILOUT)
        {
            break;
        }
    }

    float i_wormsdf = min((length(z1) * log(length(z1)) / length(dz1)) / i_XYSCALE, (length(z2) * log(length(z2)) / length(dz2)) / i_XYSCALE);

    float i_r1 = sdf_ellipsoid(pos - vec3(0, -0.07, 0), vec3(0.12, 0.15, 0.4));
    float i_r2 = sdf_ellipsoid(pos - vec3(0, -0.111, -0.24), vec3(0.18, 0.16, 0.06));
    float i_r3 = sdf_ellipsoid(pos - vec3(0, -0.103, -0.08), vec3(0.2, 0.18, 0.06));
    float i_r4 = sdf_ellipsoid(pos - vec3(0, -0.103, 0.08), vec3(0.2, 0.18, 0.06));
    float i_r5 = sdf_ellipsoid(pos - vec3(0, -0.111, 0.24), vec3(0.18, 0.16, 0.06));
    float i_innercarvesdf = min(min(min(i_r1, i_r2), min(i_r3, i_r4)), i_r5);
    float i_spheresdf = length(pos + vec3(0, -0.2, 0)) - 0.62;

    float i_endspheresdf = length(pos + vec3(0, 0.14, 0.63)) - 0.12;
    float i_startspheresdf = length(pos + vec3(0, 0.18, -0.63)) - 0.1;
    float i_mouthandass = min(i_startspheresdf, i_endspheresdf);
    return max(max(i_wormsdf, -i_innercarvesdf), -i_mouthandass);
}

float sdf(vec3 pos, int idx)
{
    // Regular box SDF.
    vec3 rpos = pos - g_world[idx * 2];
    vec3 diff = abs(rpos) - g_world[idx * 2 + 1];
    float ret = length(max(diff, 0)) + min(max(diff.x, max(diff.y, diff.z)), 0);

    if(idx == i_world_element_index_tardigrade)
    {
        float i_tolerance_tardigrade = 0.01;
        if(ret < i_tolerance_tardigrade)
        {
            return sdf_tardigrade(rpos);
        }
        return ret + i_tolerance_tardigrade;
    }

    vec2 grid = abs(max(fract(rpos.xz * 3.5) - 0.95, 0) - 0.025);
    float tiling = (grid.x + grid.y) * 0.4;

    if(idx == 2)
    {
        // Inline box frame.
        vec3 q = abs(diff + i_frame_width) - i_frame_width;
        ret = min(min(
                    length(max(vec3(diff.x, q.y, q.z), 0)) + min(max(diff.x, max(q.y, q.z)), 0),
                    length(max(vec3(q.x, diff.y, q.z), 0)) + min(max(q.x, max(diff.y, q.z)), 0)),
                length(max(vec3(q.x, q.y, diff.z), 0)) + min(max(q.x, max(q.y, diff.z)), 0));
    }
    if(idx == 0)
    {
        tiling *= 0.7 + smoothstep(0.95, 1.0, abs(rpos.y)) * 0.3;
        ret *= -1;
    }
    return ret - tiling;
}

float sdf(vec3 pos, float count, out int hit_idx)
{
    float ret = 9;

    for(int ii = 0; ii < int(count); ++ii)
    {
        if(((1 << ii) & g_object_mask) != 0)
        {
            float sdist = sdf(pos, ii);
            if(sdist < ret)
            {
                hit_idx = ii;
                ret = sdist;
            }
        }
    }

    return ret;
}

float raycast(vec3 origin, vec3 forward, vec3 settings, out vec3 hit, out vec3 nor, out int hit_idx)
{
    vec3 diff = vec3(0.01, 0, 0);
    float dist = 0;
    float ray = 0;
    float min_before_hit = 9;

    for(int ii = 0; ii < 66; ++ii)
    {
        hit = origin + forward * dist;
        dist = sdf(hit, settings.x, hit_idx);
        if(dist < 0)
        {
            break;
        }
        min_before_hit = min(min_before_hit, dist);
        origin = hit;
        dist = dist * (1 + settings.y) + settings.y;
        ray += dist;
        if(ray > settings.z)
        {
            hit_idx = i_world_elements_concrete;
            break;
        }
    }

    if(dist < 0)
    {
        for(int ii = 0; ii < 11; ++ii)
        {
            vec3 mid = (hit + origin) / 2;
            dist = sdf(mid, hit_idx);
            if(dist < 0.0)
            {
                hit = mid;
            }
            else
            {
                origin = mid;
            }
        }
    }

    nor = normalize(vec3(sdf(hit + diff.xyy, hit_idx), sdf(hit + diff.yxy, hit_idx), sdf(hit + diff.yyx, hit_idx)));

    return min_before_hit;
}

vec3 turbulence(vec3 pos, float pmul, float smul)
{
    vec3 ret = vec3(0.0);
    float scale = 1.0;
    for(int ii = 0; ii < 4; ++ii)
    {
        pos = rot * pos;
        ret += texture(u_tex, pos).xyz * scale;
        pos *= pmul;
        scale *= smul;
    }
    return ret;
}

void main()
{
    output_color = vec4(9, 9, 9, 1);

    vec3 input_position;
    vec3 input_forward;
    vec3 i_input_up = vec3(0, 1, 0);

#if defined(USE_LD)
    if(u_control[3].x == 0.0)
    {
        g_object_mask = 7 +
            ((u_control[4].x > 0.0) ? 8 : 0) +
            ((u_control[4].y > 0.0) ? 16 : 0) +
            ((u_control[4].z > 0.0) ? 32 : 0) +
            ((u_control[5].x > 0.0) ? 64 : 0) +
            ((u_control[5].y > 0.0) ? 128 : 0) +
            ((u_control[5].z > 0.0) ? 256 : 0) +
            ((u_control[6].x > 0.0) ? 512 : 0) +
            ((u_control[6].y > 0.0) ? 1024 : 0);
        rot = mat3(u_control[7].x, u_control[7].y, u_control[7].z,
                u_control[8].x, u_control[8].y, u_control[8].z,
                u_control[9].x, u_control[9].y, u_control[9].z);
        input_position = u_control[0];
        input_forward = normalize(u_control[1]);
        i_input_up = normalize(u_control[2]);
    }
    else
    {
#endif
        float ctime = u_time;
        float i_time_mul = 172 * 2048;
        for(int ii = 0; true; ii += 5)
        {
            vec3 settings = g_camera[ii];
            float ftime = settings.x * i_time_mul;
            float ratio = ctime / ftime;
            if(ratio < 1)
            {
                input_position = mix(g_camera[ii + 1], g_camera[ii + 3], ratio);
                input_forward = normalize(mix(g_camera[ii + 2], g_camera[ii + 4], ratio));
                g_world[4].z = g_world[18].z = mix(settings.y, settings.z, ratio);
                break;
            }
            ctime -= ftime;
        }
        ctime = u_time;
        for(int ii = 0; true; ii += 5)
        {
            float ftime = g_masking[ii] * i_time_mul;
            float ratio = ctime / ftime;
            if(ratio < 1)
            {
                int ormask = int(round(mix(g_masking[ii + 2], g_masking[ii + 3], ratio)));
                int shift = g_masking[ii + 4];
                int andmask = int(mix(u_time & 511, 511, (shift > 0) ? ratio : 1 - ratio));
                g_object_mask |= g_masking[ii + 1] | ((andmask & ormask) << abs(shift));
                break;
            }
            ctime -= ftime;
        }
#if defined(USE_LD)
    }
#endif

#if defined(USE_LD)
    vec2 aspect = g_position;
    if(u_control[3].y > u_control[3].z)
    {
        aspect.x *= u_control[3].y / u_control[3].z;
    }
    else
    {
        aspect.y /= u_control[3].y / u_control[3].z;
    }
    aspect *= 0.7;
#else
    vec2 aspect = g_position * vec2(1.24, 0.7); // 16:9 * 0.7
#endif

    vec3 right = cross(input_forward, i_input_up);
    input_forward = normalize(aspect.x * right + aspect.y * normalize(cross(right, input_forward)) + input_forward);

    // First hit something real.
    vec4 blend = vec4(0);
    vec3 hit;
    vec3 nor;
    float i_hit_step_forward = 0.003;
    int mask_backup = g_object_mask;
    int hit_idx;
    raycast(input_position, input_forward, vec3(i_world_elements_total, i_hit_step_forward, 9), hit, nor, hit_idx);

#if 1
    // Re-test if hit the scan but disabling it.
    while(hit_idx >= i_world_element_index_scan)
    {
        g_object_mask ^= 1 << hit_idx;

        float i_str = sdf_tardigrade(hit - g_world[20]);
        float intensity = min(max(i_str, 0) * 15, 0.5);
        vec3 i_turb1 = turbulence(hit.zxy * u_time * 0.00007, 0.5, 2.1);
        vec3 i_turb2 = turbulence(hit * intensity * 1.6, 0.3, 1.6) / 99;
        vec4 mixer = (hit_idx == i_world_element_index_scan) ?
            vec4(vec3(5, 5, 7) * sqrt(abs(i_turb1 + i_turb2)), 1) * intensity :
            vec4(9, 9, 22, 0.9) * smoothstep(-0.9, 0.0, dot(input_forward, nor));

        blend = blend * (1 - mixer.a) + mixer;

        raycast(input_position, input_forward, vec3(i_world_elements_total, i_hit_step_forward, 9), hit, nor, hit_idx);
    }
#else
    if(hit_idx == i_world_element_index_tardigrade)
    {
        hit_idx = 1;
    }
#endif

    // Calculate total lighting to hit here.
    if(hit_idx < i_world_elements_concrete)
    {
        // Disable the object we just hit.
        g_object_mask = mask_backup ^ 1 << hit_idx;

        vec3 rpos = hit - g_world[hit_idx * 2];
        float mixer = smoothstep(1.96, 1.97, max(abs(rpos.x), abs(rpos.z)));
        float rsqr = mix(0.2, 1.0, mixer);

        nor = mix(normalize(turbulence(rpos * vec3(7), 0.4, 1.2) + nor * 55),
                normalize(turbulence(rpos * vec3(5, 0.2, 5), 0.5, 1.7) + turbulence(rpos * vec3(13, 11, 13), 0.7, 2.1) + nor * 55),
                mixer);
        vec3 col = mix(vec3(0.8, 0.9, 1),
                vec3(0.95, 0.95, 1),
                mixer);
        vec3 lit = vec3(0);

        for(int ii = i_world_elements_concrete; ii < i_world_element_index_tardigrade; ++ii)
        {
            if(((1 << ii) & g_object_mask) != 0)
            {
                vec3 diff = g_world[ii * 2] - hit;
                float dist = length(diff);
                vec3 dir = diff / dist;

                // Reuse 'diff', 'right' and 'hit_idx' as they are no longer used.
                float min_before = raycast(hit, dir, vec3(i_world_elements_concrete, 0.02, dist), right, diff, hit_idx);

                // BRDF time.
                if(hit_idx >= i_world_elements_concrete)
                {
                    // Settings (inline)
                    float i_diffuse_reflectance = 1.0;
                    float i_fresnel = 0.9;

                    // Common stuff (vec3).
                    vec3 hh = normalize(dir - input_forward);
                    vec3 i_light_color = (ii == i_world_element_index_scan) ? vec3(0.25, 0.25, 0.33) : vec3(0.25);

                    // Common stuff (float).
                    float ndotl = max(dot(nor, dir), 0);
                    float ndotv = max(dot(nor, -input_forward), 0);
                    float ndoth = max(dot(nor, hh), 0);
                    float i_ldoth = max(dot(dir, hh), 0);

                    // Diffuse part.
                    float anglel = acos(ndotl);
                    float anglev = acos(ndotv);

                    float i_A = 1 - rsqr / (rsqr + 0.3) / 2; // addition constant should be 0.33
                    float i_B = rsqr / (rsqr + 0.09) / 2; // multiplier constant should be 0.45
                    float i_sina = sin(max(anglel, anglev));
                    float i_tanb = tan(min(anglel, anglev));
                    vec3 i_projected_ll = normalize(dir - ndotl * nor);
                    vec3 i_projected_vv = normalize(-input_forward - ndotv * nor);
                    float i_cos_phi_diff = max(dot(i_projected_ll, i_projected_vv), 0);

                    float i_diffuse = i_diffuse_reflectance * ndotl * (i_A + i_B * i_cos_phi_diff * i_sina * i_tanb);

                    // Specular part.
                    float alphaSqr = rsqr * rsqr;
                    float k = rsqr / 2;
                    float denom = ndoth * ndoth * (alphaSqr - 1) + 1;
                    float i_D = alphaSqr / (g_pi * denom * denom);
                    float i_F = i_fresnel + (1 - i_fresnel) * pow(1 - i_ldoth, 5);

                    float i_specular = ndotl * i_D * i_F / ((1 - ndotl * k + k) * (1 - ndotv * k + k));

                    lit += (i_diffuse + i_specular) * (i_light_color * smoothstep(0.0, 0.3, min_before * 3));
                }
            }
        }
        output_color.xyz = mix(lit * col, blend.xyz, blend.w);
    }
}
