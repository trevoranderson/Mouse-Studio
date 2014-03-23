Mouse-Studio
===========
Project Status:
You are free to use/fork/pull-request. However, I still consider this project to be in early stages,
so function names and specification are likely to change.

To use:
Put the Mat Point and MouseMovement class files into your working directory. Instantiate a MouseMovement objecct.
The stuff generally does what it says on the box.

Spec: (Often out of date but w/e)

Construction:
You can create an empty MouseMovement or instantiate one with a file.

IMPORTANT: File Loading optimistically assumes the file is a previously saved mousemovement, and will
interperet the file as a mousemovement, even if it isn't.

Recording:
Specify an amount of time to record, possibly with a resolution in points per second.
The resolution only stipulates the amount of storage used by the recording. The recording function will
partially fill a vector based on the time you specify, and use linear interpolation to fill the gaps in data.

Play:
Runs similarly to Recording, but uses SetCursorPos instead of get. Play functions do not change your underlying data.
However, if necessary to keep resolution up, they will resize a copy of the data if you provide a time that is longer
than the original, or play it over a larger area.

Load/Save:
Provides an interface for saving a recorded mousemovement, or loading from files. The files can have any extension,
it only cares about the byte data. 