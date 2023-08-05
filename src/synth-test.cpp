/// \file
/// Synth test code.

//######################################
// Include #############################
//######################################

#include "dnload.h"

#if !defined(USE_LD)
#include "synth.h"
#endif

//######################################
// Define ##############################
//######################################

/// Size of one sample in bytes.
#define AUDIO_SAMPLE_SIZE 4

/// \cond
#if (4 == AUDIO_SAMPLE_SIZE)
#define AUDIO_SAMPLE_TYPE_SDL AUDIO_F32SYS
typedef float sample_t;
#elif (2 == AUDIO_SAMPLE_SIZE)
#define AUDIO_SAMPLE_TYPE_SDL AUDIO_S16SYS
typedef int16_t sample_t;
#elif (1 == AUDIO_SAMPLE_SIZE)
#define AUDIO_SAMPLE_TYPE_SDL AUDIO_U8
typedef uint8_t sample_t;
#else
#error "invalid audio sample size"
#endif
#define AUDIO_POSITION_SHIFT (9 - (4 / sizeof(sample_t)))
/// \endcond

/// Audio channels.
#define AUDIO_CHANNELS 2

/// Audio samplerate.
#define AUDIO_SAMPLERATE 44100

/// Audio byterate.
#define AUDIO_BYTERATE (AUDIO_CHANNELS * AUDIO_SAMPLERATE * AUDIO_SAMPLE_SIZE)

/// Intro length (in bytes of audio).
#define INTRO_LENGTH (96 * AUDIO_BYTERATE)

/// Output file.
#define SYNTH_TEST_OUTPUT_FILE "tomograph.wav"

/// Debug (bytebeat) mode on/off.
#define SYNTH_DEBUG 0

//######################################
// Global data #########################
//######################################

/// "Enough". Could add exact definition in generated synth.h, but this gets us started
#if !defined(SYNTH_WORK_SIZE)
#define SYNTH_WORK_SIZE (16 * 1024 * 1024)
#endif

/// Audio buffer length (should be larger than intro length for safety).
#define AUDIO_BUFFER_LENGTH (INTRO_LENGTH * 9 / 8)

/// Audio buffer for output and work.
static uint8_t g_audio_buffer[SYNTH_WORK_SIZE + AUDIO_BUFFER_LENGTH];

/// Audio output buffer.
static uint8_t* g_audio_buffer_song = g_audio_buffer + SYNTH_WORK_SIZE;

//######################################
// Functions ###########################
//######################################

/// Synth wrapper.
/// \param output Output buffer to write to.
/// \param samples Number of bytes to write.
static void synth_func(uint8_t* output, unsigned bytes)
{
#if defined(USE_LD) || (defined(SYNTH_DEBUG) && SYNTH_DEBUG)
    {
        unsigned ii;

        // Example by "bst", taken from "Music from very short programs - the 3rd iteration" by viznut.
        for(ii = 0; ((bytes / sizeof(sample_t) / AUDIO_CHANNELS / 10) > ii); ++ii)
        {
            uint8_t sample = static_cast<uint8_t>(
                    static_cast<int>(ii / 70000000 * ii * ii + ii) % 127 |
                    ii >> 4 | ii >> 5 | (ii % 127 + (ii >> 17)) | ii
                    );
            float flt_sample = (static_cast<float>(sample) / 127.5f) - 1.0f;
            float* flt_audio_buffer = reinterpret_cast<float*>(g_audio_buffer_song);
            flt_audio_buffer[(ii * 10) + 0] = flt_sample;
            flt_audio_buffer[(ii * 10) + 1] = flt_sample;
            flt_audio_buffer[(ii * 10) + 2] = flt_sample;
            flt_audio_buffer[(ii * 10) + 3] = flt_sample;
            flt_audio_buffer[(ii * 10) + 4] = flt_sample;
            flt_audio_buffer[(ii * 10) + 5] = flt_sample;
            flt_audio_buffer[(ii * 10) + 6] = flt_sample;
            flt_audio_buffer[(ii * 10) + 7] = flt_sample;
            flt_audio_buffer[(ii * 10) + 8] = flt_sample;
            flt_audio_buffer[(ii * 10) + 9] = flt_sample;
        }
    }
#else
  synth(reinterpret_cast<float*>(output), bytes / sizeof(sample_t));
#endif

  // Write raw .wav file if necessary.
#if defined(SYNTH_TEST_OUTPUT_FILE)
  {
    SF_INFO info;
    info.samplerate = AUDIO_SAMPLERATE;
    info.channels = AUDIO_CHANNELS;
    info.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

    SNDFILE *outfile = dnload_sf_open(SYNTH_TEST_OUTPUT_FILE, SFM_WRITE, &info);
    sf_count_t write_count = INTRO_LENGTH / AUDIO_CHANNELS / AUDIO_SAMPLE_SIZE;
    dnload_sf_writef_float(outfile, reinterpret_cast<float*>(output), write_count);
    dnload_sf_close(outfile);
  }
#endif
}

//######################################
// main / _start #######################
//######################################

#if defined(USE_LD)
/// \brief Intro body function.
///
/// \param screen_w Screen width.
/// \param screen_h Screen height.
/// \param flag_fullscreen Fullscreen toggle.
/// \param flag_record Record toggle.
void main(int argc, char** argv)
#else
void _start()
#endif
{
    dnload();

    dnload_puts("synth_func()");
    synth_func(g_audio_buffer_song, INTRO_LENGTH / sizeof(sample_t));

#if !defined(USE_LD)
    dnload_puts("asm_exit()");
    asm_exit();
#endif
}

//######################################
// End #################################
//######################################

