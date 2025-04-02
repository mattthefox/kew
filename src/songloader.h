#include <glib.h>
#include <gio/gio.h>
#include <sys/wait.h>
#include <unistd.h>
#include "appstate.h"
#include "tagLibWrapper.h"
#include "cache.h"
#include "imgfunc.h"
#include "file.h"
#include "sound.h"
#include "soundcommon.h"
#include "utils.h"

#ifndef MAXPATHLEN
#define MAXPATHLEN 4096
#endif

#ifndef KEYVALUEPAIR_STRUCT
#define KEYVALUEPAIR_STRUCT

typedef struct
{
        char *key;
        char *value;
} KeyValuePair;

#endif

#ifndef TAGSETTINGS_STRUCT
#define TAGSETTINGS_STRUCT

#define METADATA_MAX_LENGTH 256

        typedef struct
        {
                char title[METADATA_MAX_LENGTH];
                char artist[METADATA_MAX_LENGTH];
                char album_artist[METADATA_MAX_LENGTH];
                char album[METADATA_MAX_LENGTH];
                char date[METADATA_MAX_LENGTH];
                double replaygainTrack;
                double replaygainAlbum;
        } TagSettings;

#endif

#ifndef SONGDATA_STRUCT
#define SONGDATA_STRUCT

typedef struct
{
        gchar *trackId;
        char filePath[MAXPATHLEN];
        char coverArtPath[MAXPATHLEN];
        unsigned char red;
        unsigned char green;
        unsigned char blue;
        unsigned char red2;
        unsigned char green2;
        unsigned char blue2;
        TagSettings *metadata;
        unsigned char *cover;
        int coverWidth;
        int coverHeight;
        double duration;
        bool hasErrors;
} SongData;

#endif

SongData *loadSongData(char *filePath, AppState *state);

void unloadSongData(SongData **songdata, AppState *state);
