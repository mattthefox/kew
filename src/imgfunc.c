#include "imgfunc.h"
#include "term.h"

/*

imgfunc.c

 Related to displaying an image in the terminal.

*/

// Disable some warnings for stb headers.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wstrict-overflow"
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include <stb_image_resize2.h>
#pragma GCC diagnostic pop

/*

chafafunc.c

 Functions related to printing images to the terminal with chafa.

*/

/* Include after chafa.h for G_OS_WIN32 */
#ifdef G_OS_WIN32
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#endif
#include <io.h>
#else
#include <sys/ioctl.h> /* ioctl */
#endif

#define MACRO_STRLEN(s) (sizeof(s) / sizeof(s[0]))
#define MAX_COLOR_COUNT 10  // Max number of color clusters to consider

typedef struct
{
        gint width_cells, height_cells;
        gint width_pixels, height_pixels;
} TermSize;

char scale[] = "$@&B%8WM#ZO0QoahkbdpqwmLCJUYXIjft/\\|()1{}[]l?zcvunxr!<>i;:*-+~_,\"^`'. ";
unsigned int brightness_levels = MACRO_STRLEN(scale) - 2;

static void detect_terminal(ChafaTermInfo **term_info_out, ChafaCanvasMode *mode_out, ChafaPixelMode *pixel_mode_out)
{
        ChafaCanvasMode mode;
        ChafaPixelMode pixel_mode;
        ChafaTermInfo *term_info;
        gchar **envp;

        /* Examine the environment variables and guess what the terminal can do */

        envp = g_get_environ();
        term_info = chafa_term_db_detect(chafa_term_db_get_default(), envp);

        /* See which control sequences were defined, and use that to pick the most
         * high-quality rendering possible */

        if (chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_BEGIN_KITTY_IMMEDIATE_IMAGE_V1))
        {
                pixel_mode = CHAFA_PIXEL_MODE_KITTY;
                mode = CHAFA_CANVAS_MODE_TRUECOLOR;
        }
        else if (chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_BEGIN_SIXELS))
        {
                pixel_mode = CHAFA_PIXEL_MODE_SIXELS;
                mode = CHAFA_CANVAS_MODE_TRUECOLOR;
        }
        else
        {
                pixel_mode = CHAFA_PIXEL_MODE_SYMBOLS;

                if (chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_FGBG_DIRECT) && chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_FG_DIRECT) && chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_BG_DIRECT))
                        mode = CHAFA_CANVAS_MODE_TRUECOLOR;
                else if (chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_FGBG_256) && chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_FG_256) && chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_BG_256))
                        mode = CHAFA_CANVAS_MODE_INDEXED_240;
                else if (chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_FGBG_16) && chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_FG_16) && chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_SET_COLOR_BG_16))
                        mode = CHAFA_CANVAS_MODE_INDEXED_16;
                else if (chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_INVERT_COLORS) && chafa_term_info_have_seq(term_info, CHAFA_TERM_SEQ_RESET_ATTRIBUTES))
                        mode = CHAFA_CANVAS_MODE_FGBG_BGFG;
                else
                        mode = CHAFA_CANVAS_MODE_FGBG;
        }

        /* Hand over the information to caller */

        *term_info_out = term_info;
        *mode_out = mode;
        *pixel_mode_out = pixel_mode;

        /* Cleanup */

        g_strfreev(envp);
}

static void
get_tty_size(TermSize *term_size_out)
{
        TermSize term_size;

        term_size.width_cells = term_size.height_cells = term_size.width_pixels = term_size.height_pixels = -1;

#ifdef G_OS_WIN32
        {
                HANDLE chd = GetStdHandle(STD_OUTPUT_HANDLE);
                CONSOLE_SCREEN_BUFFER_INFO csb_info;

                if (chd != INVALID_HANDLE_VALUE && GetConsoleScreenBufferInfo(chd, &csb_info))
                {
                        term_size.width_cells = csb_info.srWindow.Right - csb_info.srWindow.Left + 1;
                        term_size.height_cells = csb_info.srWindow.Bottom - csb_info.srWindow.Top + 1;
                }
        }
#else
        {
                struct winsize w;
                gboolean have_winsz = FALSE;

                if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) >= 0 || ioctl(STDERR_FILENO, TIOCGWINSZ, &w) >= 0 || ioctl(STDIN_FILENO, TIOCGWINSZ, &w) >= 0)
                        have_winsz = TRUE;

                if (have_winsz)
                {
                        term_size.width_cells = w.ws_col;
                        term_size.height_cells = w.ws_row;
                        term_size.width_pixels = w.ws_xpixel;
                        term_size.height_pixels = w.ws_ypixel;
                }
        }
#endif

        if (term_size.width_cells <= 0)
                term_size.width_cells = -1;
        if (term_size.height_cells <= 2)
                term_size.height_cells = -1;

        /* If .ws_xpixel and .ws_ypixel are filled out, we can calculate
         * aspect information for the font used. Sixel-capable terminals
         * like mlterm set these fields, but most others do not. */

        if (term_size.width_pixels <= 0 || term_size.height_pixels <= 0)
        {
                term_size.width_pixels = -1;
                term_size.height_pixels = -1;
        }

        *term_size_out = term_size;
}

static void
tty_init(void)
{
#ifdef G_OS_WIN32
        {
                HANDLE chd = GetStdHandle(STD_OUTPUT_HANDLE);

                saved_console_output_cp = GetConsoleOutputCP();
                saved_console_input_cp = GetConsoleCP();

                /* Enable ANSI escape sequence parsing etc. on MS Windows command prompt */

                if (chd != INVALID_HANDLE_VALUE)
                {
                        if (!SetConsoleMode(chd,
                                            ENABLE_PROCESSED_OUTPUT | ENABLE_WRAP_AT_EOL_OUTPUT | ENABLE_VIRTUAL_TERMINAL_PROCESSING | DISABLE_NEWLINE_AUTO_RETURN))
                                win32_stdout_is_file = TRUE;
                }

                /* Set UTF-8 code page I/O */

                SetConsoleOutputCP(65001);
                SetConsoleCP(65001);
        }
#endif
}

static GString *
convert_image(const void *pixels, gint pix_width, gint pix_height,
              gint pix_rowstride, ChafaPixelType pixel_type,
              gint width_cells, gint height_cells,
              gint cell_width, gint cell_height)
{
        ChafaTermInfo *term_info;
        ChafaCanvasMode mode;
        ChafaPixelMode pixel_mode;
        ChafaSymbolMap *symbol_map;
        ChafaCanvasConfig *config;
        ChafaCanvas *canvas;
        GString *printable;

        detect_terminal(&term_info, &mode, &pixel_mode);

        /* Specify the symbols we want */

        symbol_map = chafa_symbol_map_new();
        chafa_symbol_map_add_by_tags(symbol_map, CHAFA_SYMBOL_TAG_BLOCK);

        /* Set up a configuration with the symbols and the canvas size in characters */

        config = chafa_canvas_config_new();
        chafa_canvas_config_set_canvas_mode(config, mode);
        chafa_canvas_config_set_pixel_mode(config, pixel_mode);
        chafa_canvas_config_set_geometry(config, width_cells, height_cells);
        chafa_canvas_config_set_symbol_map(config, symbol_map);

        if (cell_width > 0 && cell_height > 0)
        {
                /* We know the pixel dimensions of each cell. Store it in the config. */

                chafa_canvas_config_set_cell_geometry(config, cell_width, cell_height);
        }

        /* Create canvas */

        canvas = chafa_canvas_new(config);

        /* Draw pixels to the canvas */

        chafa_canvas_draw_all_pixels(canvas,
                                     pixel_type,
                                     pixels,
                                     pix_width,
                                     pix_height,
                                     pix_rowstride);

        /* Build printable string */
        printable = chafa_canvas_print(canvas, term_info);

        /* Clean up and return */

        chafa_canvas_unref(canvas);
        chafa_canvas_config_unref(config);
        chafa_symbol_map_unref(symbol_map);
        chafa_term_info_unref(term_info);
        canvas = NULL;
        config = NULL;
        symbol_map = NULL;
        term_info = NULL;
        return printable;
}

unsigned char *getBitmap(const char *image_path, int *width, int *height)
{
    if (image_path == NULL)
        return NULL;

    int channels;

    // Load the image with original dimensions first
    unsigned char *image = stbi_load(image_path, width, height, &channels, 4); // Force 4 channels (RGBA)
    if (!image)
    {
        fprintf(stderr, "Failed to load image: %s\n", image_path);
        return NULL;
    }

    // Ensure width and height are equal by cropping from the center
    int min_dim = (*width < *height) ? *width : *height;

    // Calculate the cropping area
    int crop_x = (*width - min_dim) / 2;
    int crop_y = (*height - min_dim) / 2;

    // Create a new array to store the square image
    unsigned char *square_image = (unsigned char *)malloc(min_dim * min_dim * 4);  // RGBA, 4 bytes per pixel
    if (!square_image)
    {
        fprintf(stderr, "Failed to allocate memory for the square image.\n");
        stbi_image_free(image);
        return NULL;
    }

    // Copy pixels to the new square image (cropping the center)
    for (int y = 0; y < min_dim; y++)
    {
        for (int x = 0; x < min_dim; x++)
        {
            int src_index = ((crop_y + y) * (*width) + (crop_x + x)) * 4;
            int dest_index = (y * min_dim + x) * 4;
            square_image[dest_index + 0] = image[src_index + 0]; // Red
            square_image[dest_index + 1] = image[src_index + 1]; // Green
            square_image[dest_index + 2] = image[src_index + 2]; // Blue
            square_image[dest_index + 3] = image[src_index + 3]; // Alpha
        }
    }

    // Update the width and height to the square dimensions
    *width = *height = min_dim;

    // Free the original image since it's no longer needed
    stbi_image_free(image);

    return square_image;
}


float calcAspectRatio(void)
{
        TermSize term_size;
        gint cell_width = -1, cell_height = -1;

        tty_init();
        get_tty_size(&term_size);

        if (term_size.width_cells > 0 && term_size.height_cells > 0 && term_size.width_pixels > 0 && term_size.height_pixels > 0)
        {
                cell_width = term_size.width_pixels / term_size.width_cells;
                cell_height = term_size.height_pixels / term_size.height_cells;
        }

        // Set default for some terminals
        if (cell_width == -1 && cell_height == -1)
        {
                cell_width = 8;
                cell_height = 16;
        }

        return (float)cell_height / (float)cell_width;
}

void printSquareBitmapCentered(unsigned char *pixels, int width, int height, int baseHeight)
{
        if (pixels == NULL)
        {
                printf("Error: Invalid pixel data.\n");
                return;
        }

        // Use the provided width and height
        int pix_width = width;
        int pix_height = height;
        int n_channels = 4; // Assuming RGBA format

        // Validate the image dimensions
        if (pix_width == 0 || pix_height == 0)
        {
                printf("Error: Invalid image dimensions.\n");
                return;
        }

        TermSize term_size;
        GString *printable;
        gint cell_width = -1, cell_height = -1;

        tty_init();
        get_tty_size(&term_size);

        if (term_size.width_cells > 0 && term_size.height_cells > 0 &&
            term_size.width_pixels > 0 && term_size.height_pixels > 0)
        {
                cell_width = term_size.width_pixels / term_size.width_cells;
                cell_height = term_size.height_pixels / term_size.height_cells;
        }

        // Set default cell size for some terminals
        if (cell_width == -1 || cell_height == -1)
        {
                cell_width = 8;
                cell_height = 16;
        }

        // Calculate corrected width based on aspect ratio correction
        float aspect_ratio_correction = (float)cell_height / (float)cell_width;
        int correctedWidth = (int)(baseHeight * aspect_ratio_correction);

        // Convert image to a printable string using Chafa
        printable = convert_image(
            pixels,
            pix_width,
            pix_height,
            pix_width * n_channels,         // Row stride
            CHAFA_PIXEL_RGBA8_UNASSOCIATED, // Correct pixel format
            correctedWidth,
            baseHeight,
            cell_width,
            cell_height);

        // Ensure the string is null-terminated
        g_string_append_c(printable, '\0');

        // Split the printable string into lines
        const gchar *delimiters = "\n";
        gchar **lines = g_strsplit(printable->str, delimiters, -1);

        // Calculate indentation to center the image
        int indentation = ((term_size.width_cells - correctedWidth) / 2);

        // Print each line with indentation
        for (int i = 0; lines[i] != NULL; i++)
        {
                printf("\n\033[%dC%s", indentation, lines[i]);
        }

        // Free allocated memory
        g_strfreev(lines);
        g_string_free(printable, TRUE);
}

unsigned char luminanceFromRGB(unsigned char r, unsigned char g, unsigned char b)
{
        return (unsigned char)(0.2126 * r + 0.7152 * g + 0.0722 * b);
}

void checkIfBrightPixel(unsigned char r, unsigned char g, unsigned char b, bool *found)
{
        // Calc luminace and use to find Ascii char.
        unsigned char ch = luminanceFromRGB(r, g, b);

        if (ch > 80 && !(r < g + 20 && r > g - 20 && g < b + 20 && g > b - 20) && !(r > 150 && g > 150 && b > 150))
        {
                *found = true;
        }
}

float colorDistance(unsigned char r1, unsigned char g1, unsigned char b1, unsigned char r2, unsigned char g2, unsigned char b2)
{
    return sqrtf(powf(r1 - r2, 2) + powf(g1 - g2, 2) + powf(b1 - b2, 2));
}

// Get the most frequent color in the image (simple average approach)
int getCoverColor(unsigned char *pixels, int width, int height, unsigned char *r, unsigned char *g, unsigned char *b, unsigned char *r2, unsigned char *g2, unsigned char *b2)
{
    if (pixels == NULL || width <= 0 || height <= 0)
    {
        return -1;
    }

    int channels = 4;  // RGBA format
    int numPixels = width * height;

    // To hold colors (up to MAX_COLOR_COUNT most frequent colors)
    unsigned char colorBuckets[MAX_COLOR_COUNT][3] = {0};  // RGB
    int colorCounts[MAX_COLOR_COUNT] = {0};

    // Sampling step: We choose a fraction of pixels to analyze
    int sampleStep = numPixels / MAX_COLOR_COUNT; 

    for (int i = 0; i < numPixels; i += sampleStep)
    {
        int index = i * channels;
        unsigned char red = pixels[index + 0];
        unsigned char green = pixels[index + 1];
        unsigned char blue = pixels[index + 2];

        // Find where to store this color
        //int foundBucket = -1;
        for (int j = 0; j < MAX_COLOR_COUNT; j++)
        {
            // If this bucket is empty, fill it with this color
            if (colorCounts[j] == 0)
            {
                colorBuckets[j][0] = red;
                colorBuckets[j][1] = green;
                colorBuckets[j][2] = blue;
                colorCounts[j] = 1;
                //foundBucket = j;
                break;
            }
            // If this color is similar to the bucket, increment its count
            if (colorDistance(red, green, blue, colorBuckets[j][0], colorBuckets[j][1], colorBuckets[j][2]) < 50.0f)
            {
                colorBuckets[j][0] = (colorBuckets[j][0] * colorCounts[j] + red) / (colorCounts[j] + 1);
                colorBuckets[j][1] = (colorBuckets[j][1] * colorCounts[j] + green) / (colorCounts[j] + 1);
                colorBuckets[j][2] = (colorBuckets[j][2] * colorCounts[j] + blue) / (colorCounts[j] + 1);
                colorCounts[j]++;
                //foundBucket = j;
                break;
            }
        }
    }

    // Now find two most contrasting colors
    float maxDistance = 0.0f;
    int color1Idx = 0, color2Idx = 0;

    for (int i = 0; i < MAX_COLOR_COUNT; i++)
    {
        if (colorCounts[i] == 0) continue;  // Skip empty buckets
        for (int j = i + 1; j < MAX_COLOR_COUNT; j++)
        {
            if (colorCounts[j] == 0) continue;  // Skip empty buckets
            float dist = colorDistance(colorBuckets[i][0], colorBuckets[i][1], colorBuckets[i][2], colorBuckets[j][0], colorBuckets[j][1], colorBuckets[j][2]);
            if (dist > maxDistance)
            {
                maxDistance = dist;
                color1Idx = i;
                color2Idx = j;
            }
        }
    }

    // Set the two contrasting colors
    *r = colorBuckets[color1Idx][0];
    *g = colorBuckets[color1Idx][1];
    *b = colorBuckets[color1Idx][2];

    *r2 = colorBuckets[color2Idx][0];
    *g2 = colorBuckets[color2Idx][1];
    *b2 = colorBuckets[color2Idx][2];

    return 0;
}


unsigned char calcAsciiChar(PixelData *p)
{
        unsigned char ch = luminanceFromRGB(p->r, p->g, p->b);
        int rescaled = ch * brightness_levels / 256;

        return scale[brightness_levels - rescaled];
}

int convertToAscii(const char *filepath, unsigned int height)
{
        /*
        Modified, originally by Danny Burrows:
        https://github.com/danny-burrows/img_to_txt

        MIT License

        Copyright (c) 2021 Danny Burrows

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.
        */

        TermSize term_size;
        gint cell_width = -1, cell_height = -1;

        tty_init();
        get_tty_size(&term_size);

        if (term_size.width_cells > 0 && term_size.height_cells > 0 &&
            term_size.width_pixels > 0 && term_size.height_pixels > 0)
        {
                cell_width = term_size.width_pixels / term_size.width_cells;
                cell_height = term_size.height_pixels / term_size.height_cells;
        }

        float aspect_ratio_correction = (float)cell_height / (float)cell_width;
        unsigned int correctedWidth = (int)(height * aspect_ratio_correction) - 1;


        // Calculate indentation to center the image
        int indent = ((term_size.width_cells - correctedWidth) / 2);

        int rwidth, rheight, rchannels;
        unsigned char *read_data = stbi_load(filepath, &rwidth, &rheight, &rchannels, 3);

        if (read_data == NULL)
        {
                return -1;
        }

        PixelData *data;
        if (correctedWidth != (unsigned)rwidth || height != (unsigned)rheight)
        {
                // 3 * uint8 for RGB!
                unsigned char *new_data = malloc(3 * sizeof(unsigned char) * correctedWidth * height);
                stbir_resize_uint8_srgb(
                    read_data, rwidth, rheight, 0,
                    new_data, correctedWidth, height, 0, 3);

                stbi_image_free(read_data);
                data = (PixelData *)new_data;
        }
        else
        {
                data = (PixelData *)read_data;
        }

        printf("\n");
        printf("%*s", indent, "");

        for (unsigned int d = 0; d < correctedWidth * height; d++)
        {
                if (d % correctedWidth == 0 && d != 0)
                {
                        printf("\n");
                        printf("%*s", indent, "");
                }

                PixelData *c = data + d;

                printf("\033[1;38;2;%03u;%03u;%03um%c", c->r, c->g, c->b, calcAsciiChar(c));
        }

        printf("\n");

        stbi_image_free(data);
        return 0;
}

int printInAscii(const char *pathToImgFile, int height)
{
        printf("\r");

        int ret = convertToAscii(pathToImgFile,(unsigned)height);
        if (ret == -1)
                printf("\033[0m");
        return 0;
}
