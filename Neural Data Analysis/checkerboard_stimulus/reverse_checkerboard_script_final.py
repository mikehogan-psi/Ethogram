from psychopy import visual, core, event
import os
import serial

# === Constants ===
SCREEN_WIDTH, SCREEN_HEIGHT = 1280, 720
GREY_COLOR = [0, 0, 0]
SWITCH_FREQ_HZ = 0.5
SWITCH_INTERVAL = 1.0 / SWITCH_FREQ_HZ  # 2 seconds per full cycle
STIM_DURATION_MIN = 5
TOTAL_DURATION_SEC = STIM_DURATION_MIN * 60
FRAMES_PER_SECOND = 15

# === Serial Trigger Settings ===
use_openephys_serial = True # True

openephys_ser = serial.Serial('COM9', 9600, timeout=1) if use_openephys_serial else None

def send_triggers():
    try:
        if openephys_ser:
            openephys_ser.write(b't')
    except:
        pass

# === Checkerboard Setup ===
checkerboard_folder = r"C:\PhD 1st Year\Fear Conditioning Behavioural Pilot Videos and Code\checkerboard_stimulus\checkerboard_images\big"
checkerboard_files = sorted([f for f in os.listdir(checkerboard_folder) if f.endswith('.png')])

if len(checkerboard_files) < 2:
    raise ValueError("Need at least two checkerboard images.")

checkerboard_paths = [os.path.join(checkerboard_folder, f) for f in checkerboard_files]

# === Window and Stimuli ===
win = visual.Window(
    size=(SCREEN_WIDTH, SCREEN_HEIGHT),
    fullscr=True,
    color=GREY_COLOR,
    units='pix'
)

checker_images = [
    visual.ImageStim(win, image=checkerboard_paths[0], size=(SCREEN_WIDTH, SCREEN_HEIGHT)),
    visual.ImageStim(win, image=checkerboard_paths[1], size=(SCREEN_WIDTH, SCREEN_HEIGHT))
]

def play_checkerboard_loop(images, duration_sec, switch_interval):
    core.wait(0.5)  # Allow PsychoPy to settle
    win.flip()      # Sync to refresh
    frame_duration = 1.0 / FRAMES_PER_SECOND
    next_flip_time = core.getTime() + frame_duration

    total_cycles = int(duration_sec / switch_interval)
    index = 0


    # === Stimulus Alternation ===
    for _ in range(total_cycles):
        core.wait(max(0, next_flip_time - core.getTime()))
        images[index].draw()
        win.flip()
        send_triggers()
        if 'q' in event.getKeys():
            return True
        index = (index + 1) % 2
        next_flip_time += switch_interval


    return False

# === Run Stimulus Presentation ===
quit_early = play_checkerboard_loop(checker_images, TOTAL_DURATION_SEC, SWITCH_INTERVAL)

if camera_ser:
    camera_ser.close()
if openephys_ser:
    openephys_ser.close()

win.close()
core.quit()
