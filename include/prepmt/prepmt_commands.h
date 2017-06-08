#ifndef PREPMT_COMMANDS_H__
#define PREPMT_COMMANDS_H__ 1
#include "prepmt/prepmt_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif

int prepmt_commands_freePrepmtCommands(struct prepmtCommands_struct *cmds);

struct prepmtCommands_struct
    prepmt_commands_readFromIniFile(const char *iniFile,
                                     const char *section,
                                     const int nobs,
                                     const struct sacData_struct *data,
                                     int *ierr);

int prepmt_commands_modifyCommandsCharsStruct(
    const struct prepmtModifyCommands_struct options,
    const struct sacData_struct data,
    struct prepmtCommandsChars_struct *cmds);

char **prepmt_commands_modifyCommands(
    const int ncmds, const char **cmds,
    const struct prepmtModifyCommands_struct options,
    const struct sacData_struct data,
    int *ierr);
#ifdef __cplusplus
}
#endif
#endif
