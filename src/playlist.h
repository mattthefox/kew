#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE
#endif
#ifndef __USE_XOPEN_EXTENDED
#define __USE_XOPEN_EXTENDED
#endif

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file.h"

#ifndef PLAYLIST_STRUCT
#define PLAYLIST_STRUCT

typedef struct
{
        char *filePath;
        double duration;
} SongInfo;

typedef struct Node
{
        int id;
        SongInfo song;
        struct Node *next;
        struct Node *prev;
} Node;

typedef struct
{
        Node *head;
        Node *tail;
        int count; 
        pthread_mutex_t mutex;
} PlayList;

extern Node *currentSong;

#endif

extern PlayList playlist;

extern PlayList *specialPlaylist;
extern PlayList *originalPlaylist;
extern int nodeIdCounter;

Node *getListNext(Node *node);

Node *getListPrev(Node *node);

void createNode(Node **node, const char *directoryPath, int id);

void addToList(PlayList *list, Node *newNode);

Node *deleteFromList(PlayList *list, Node *node);

void deletePlaylist(PlayList *playlist);

void shufflePlaylist(PlayList *playlist);

void shufflePlaylistStartingFromSong(PlayList *playlist, Node *song);

int makePlaylist(int argc, char *argv[], bool exactSearch, const char *path);

void writeM3UFile(const char *filename, PlayList *playlist);

void loadMainPlaylist(const char *directory);

void saveSpecialPlaylist(const char *directory, bool isPlayingMain);

void savePlaylist(const char *path);

PlayList deepCopyPlayList(PlayList *originalList);

void deepCopyPlayListOntoList(PlayList *originalList, PlayList *newList);

Node *findPathInPlaylist(char *path, PlayList *playlist);

Node *findLastPathInPlaylist(char *path, PlayList *playlist);

int findNodeInList(PlayList *list, int id, Node **foundNode);
