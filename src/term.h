#ifndef TERM_H
#define TERM_H
#ifndef __USE_POSIX
#define __USE_POSIX
#endif

#include <fcntl.h>
#include <poll.h>
#include <pwd.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/param.h>
#include <sys/select.h>
#include <termios.h>
#include <unistd.h>

#ifdef __GNU__
# define _BSD_SOURCE
#endif

extern volatile sig_atomic_t resizeFlag;

void setTextColor(int color);

void setTextColorRGB(int r, int g, int b);

void getTermSize(int *width, int *height);

int getIndentation(int terminalWidth);

void setNonblockingMode(void);

void restoreTerminalMode(void);

void setDefaultTextColor(void);

int isInputAvailable(void);

void resetConsole(void);

void saveCursorPosition(void);

void restoreCursorPosition(void);

void hideCursor(void);

void showCursor(void);

void clearRestOfScreen(void);

void enableScrolling(void);

void handleResize(int sig);

void resetResizeFlag(int sig);

void initResize(void);

void disableInputBuffering(void);

void enableInputBuffering(void);

void cursorJump(int numRows);

void cursorJumpDown(int numRows);

void clearScreen(void);

int readInputSequence(char *seq, size_t seqSize);

int isFunctionKey(const char *seq);

void convertControlNotationToAscii(char *input, char *output, size_t size);

void convertAsciiToControlNotation(char input, char *output, size_t size);

#endif
