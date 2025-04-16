#ifndef SETTINGS_H
#define SETTINGS_H

#include <pwd.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <unistd.h>
#include "appstate.h"
#include "events.h"
#include "file.h"
#include "soundcommon.h"
#include "player.h"
#include "utils.h"

#ifndef MAXPATHLEN
#define MAXPATHLEN 4096
#endif

#ifndef NUM_KEY_MAPPINGS
#define NUM_KEY_MAPPINGS 64
#endif

extern AppSettings settings;

enum EventType getMouseAction(int num);

void getConfig(AppSettings *settings, UISettings *ui);

void setConfig(AppSettings *settings, UISettings *ui);

void mapSettingsToKeys(AppSettings *settings, UISettings *ui, EventMapping *mappings);

#endif
