/********************************************************
 *
 * Structure for ext_data_fetch_event() filling.
 *
 * Do not edit - automatically generated.
 */

#ifndef __GUARD_H101_WRSOFIA_EXT_H101_WRSOFIA_H__
#define __GUARD_H101_WRSOFIA_EXT_H101_WRSOFIA_H__

#ifndef __CINT__
# include <stdint.h>
#else
/* For CINT (old version trouble with stdint.h): */
# ifndef uint32_t
typedef unsigned int uint32_t;
typedef          int  int32_t;
# endif
#endif
#ifndef EXT_STRUCT_CTRL
# define EXT_STRUCT_CTRL(x)
#endif

/********************************************************
 *
 * Plain structure (layout as ntuple/root file):
 */

typedef struct EXT_STR_h101_WRSOFIA_t
{
  /* RAW */
  uint32_t TIMESTAMP_SOFIA_ID /* [0,65535] */;
  uint32_t TIMESTAMP_SOFIA_WR_T1 /* [0,65535] */;
  uint32_t TIMESTAMP_SOFIA_WR_T2 /* [0,65535] */;
  uint32_t TIMESTAMP_SOFIA_WR_T3 /* [0,65535] */;
  uint32_t TIMESTAMP_SOFIA_WR_T4 /* [0,65535] */;

} EXT_STR_h101_WRSOFIA;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_WRSOFIA_onion_t
{
  /* RAW */
  uint32_t TIMESTAMP_SOFIA_ID;
  uint32_t TIMESTAMP_SOFIA_WR_T[4];

} EXT_STR_h101_WRSOFIA_onion;

/*******************************************************/

#define EXT_STR_h101_WRSOFIA_ITEMS_INFO(ok,si,offset,struct_t,printerr) do { \
  ok = 1; \
  /* RAW */ \
  EXT_STR_ITEM_INFO_LIM(ok,si,offset,struct_t,printerr,\
                     TIMESTAMP_SOFIA_ID,              UINT32,\
                    "TIMESTAMP_SOFIA_ID",65535); \
  EXT_STR_ITEM_INFO_LIM(ok,si,offset,struct_t,printerr,\
                     TIMESTAMP_SOFIA_WR_T1,           UINT32,\
                    "TIMESTAMP_SOFIA_WR_T1",65535); \
  EXT_STR_ITEM_INFO_LIM(ok,si,offset,struct_t,printerr,\
                     TIMESTAMP_SOFIA_WR_T2,           UINT32,\
                    "TIMESTAMP_SOFIA_WR_T2",65535); \
  EXT_STR_ITEM_INFO_LIM(ok,si,offset,struct_t,printerr,\
                     TIMESTAMP_SOFIA_WR_T3,           UINT32,\
                    "TIMESTAMP_SOFIA_WR_T3",65535); \
  EXT_STR_ITEM_INFO_LIM(ok,si,offset,struct_t,printerr,\
                     TIMESTAMP_SOFIA_WR_T4,           UINT32,\
                    "TIMESTAMP_SOFIA_WR_T4",65535); \
  \
} while (0);

/********************************************************
 *
 * For internal use by the network data reader:
 * (version checks, etc)
 */

typedef struct EXT_STR_h101_WRSOFIA_layout_t
{
  uint32_t _magic;
  uint32_t _size_info;
  uint32_t _size_struct;
  uint32_t _size_struct_onion;
  uint32_t _pack_list_items;

  uint32_t _num_items;
  struct {
    uint32_t _offset;
    uint32_t _size;
    uint32_t _xor;
    const char *_name;
  } _items[1];
  uint32_t _pack_list[7];
} EXT_STR_h101_WRSOFIA_layout;

#define EXT_STR_h101_WRSOFIA_LAYOUT_INIT { \
  0x57e65c96, \
  sizeof(EXT_STR_h101_WRSOFIA_layout), \
  sizeof(EXT_STR_h101_WRSOFIA), \
  sizeof(EXT_STR_h101_WRSOFIA_onion), \
  7, \
  1, \
  { \
    { 0, sizeof(EXT_STR_h101_WRSOFIA), 0x2b4a0aa0, "h101_WRSOFIA" }, \
  }, \
  { \
    0x40000000, 0x40000004, 0x40000008, 0x4000000c, \
    0x40000010, 0x40000014, 0x40000018, \
  } \
};

#endif/*__GUARD_H101_WRSOFIA_EXT_H101_WRSOFIA_H__*/

/*******************************************************/
