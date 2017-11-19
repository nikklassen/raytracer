#![windows_subsystem = "windows"]
extern crate gdi32;
extern crate kernel32;
extern crate user32;
extern crate winapi;
extern crate rayon;

mod tracer;
mod canvas;

use std::sync::atomic::{AtomicIsize, ATOMIC_ISIZE_INIT, Ordering};
use std::ffi::OsStr;
use std::os::windows::ffi::OsStrExt;
use std::iter::once;
use std::mem;
use std::ptr::{null_mut, null};
use std::io::Error;

use gdi32::StretchDIBits;
use winapi::{c_int, c_void, BITMAPINFO, BITMAPINFOHEADER, BI_RGB, DIB_RGB_COLORS, HDC, HWND,
             PAINTSTRUCT, SRCCOPY};
use winapi::minwindef::{LPARAM, LRESULT, UINT, WPARAM};
use winapi::winuser::{CS_HREDRAW, CS_OWNDC, CS_VREDRAW, CW_USEDEFAULT, MSG, WM_PAINT, WNDCLASSW,
                      WS_OVERLAPPEDWINDOW, WS_VISIBLE, VK_LEFT, VK_RIGHT, VK_UP, VK_DOWN, WM_KEYUP};
use user32::{BeginPaint, CreateWindowExW, DefWindowProcW, DispatchMessageW, EndPaint, GetMessageW,
             RegisterClassW, TranslateMessage, InvalidateRect};
use kernel32::GetModuleHandleW;

use canvas::Canvas;

fn win32_string(value: &str) -> Vec<u16> {
    OsStr::new(value).encode_wide().chain(once(0)).collect()
}

struct Window {
    handle: HWND,
}

const WIN_WIDTH: u32 = 600;
const WIN_HEIGHT: u32 = 600;

fn create_window(name: &str, title: &str) -> Result<Window, Error> {
    let name = win32_string(name);
    let title = win32_string(title);

    unsafe {
        let hinstance = GetModuleHandleW(null_mut());
        let wnd_class = WNDCLASSW {
            style: CS_OWNDC | CS_HREDRAW | CS_VREDRAW,
            lpfnWndProc: Some(message_handler),
            hInstance: hinstance,
            lpszClassName: name.as_ptr(),
            cbClsExtra: 0,
            cbWndExtra: 0,
            hIcon: null_mut(),
            hCursor: null_mut(),
            hbrBackground: null_mut(),
            lpszMenuName: null_mut(),
        };

        RegisterClassW(&wnd_class);

        let handle = CreateWindowExW(
            0,
            name.as_ptr(),
            title.as_ptr(),
            WS_OVERLAPPEDWINDOW | WS_VISIBLE,
            CW_USEDEFAULT,
            CW_USEDEFAULT,
            WIN_WIDTH as c_int,
            (WIN_HEIGHT + 50) as c_int,
            null_mut(),
            null_mut(),
            hinstance,
            null_mut(),
        );

        if handle.is_null() {
            Err(Error::last_os_error())
        } else {
            Ok(Window { handle })
        }
    }
}

static X_ROTATION: AtomicIsize = ATOMIC_ISIZE_INIT;
static Y_ROTATION: AtomicIsize = ATOMIC_ISIZE_INIT;

unsafe extern "system" fn message_handler(
    hwnd: HWND,
    u_msg: UINT,
    w_param: WPARAM,
    l_param: LPARAM,
) -> LRESULT {
    if u_msg == WM_PAINT {
        let mut canvas = Canvas::new(WIN_WIDTH as i32, WIN_HEIGHT as i32);
        tracer::draw(&mut canvas, X_ROTATION.load(Ordering::SeqCst) as f32, Y_ROTATION.load(Ordering::SeqCst) as f32);

        let mut info: BITMAPINFO = mem::uninitialized();
        info.bmiHeader.biBitCount = 24;
        info.bmiHeader.biWidth = WIN_WIDTH as i32;
        info.bmiHeader.biHeight = WIN_HEIGHT as i32;
        info.bmiHeader.biPlanes = 1;
        info.bmiHeader.biSize = mem::size_of::<BITMAPINFOHEADER>() as u32;
        info.bmiHeader.biSizeImage = 0;
        info.bmiHeader.biCompression = BI_RGB;

        let mut ps: PAINTSTRUCT = mem::uninitialized();
        let hdc: HDC = BeginPaint(hwnd, &mut ps as *mut PAINTSTRUCT);
        StretchDIBits(
            hdc,
            0,
            0,
            WIN_WIDTH as c_int,
            WIN_HEIGHT as c_int,
            0,
            0,
            WIN_WIDTH as c_int,
            WIN_HEIGHT as c_int,
            canvas.buffer.as_slice().as_ptr() as *const c_void,
            &info,
            DIB_RGB_COLORS,
            SRCCOPY,
        );
        EndPaint(hwnd, &ps as *const PAINTSTRUCT);
        0
    } else {
        DefWindowProcW(hwnd, u_msg, w_param, l_param)
    }
}

fn handle_message(window: &mut Window) -> bool {
    unsafe {
        let mut message: MSG = mem::uninitialized();
        if GetMessageW(&mut message as *mut MSG, window.handle, 0, 0) > 0 {
            if message.message == WM_KEYUP {
                if message.wParam == (VK_LEFT as u64) {
                    Y_ROTATION.fetch_sub(1, Ordering::SeqCst);
                } else if message.wParam == (VK_RIGHT as u64) {
                    Y_ROTATION.fetch_add(1, Ordering::SeqCst);
                } else if message.wParam == (VK_UP as u64) {
                    X_ROTATION.fetch_sub(1, Ordering::SeqCst);
                } else if message.wParam == (VK_DOWN as u64) {
                    X_ROTATION.fetch_add(1, Ordering::SeqCst);
                }
                InvalidateRect(message.hwnd, null(), 1);
            } else {
                TranslateMessage(&message as *const MSG);
                DispatchMessageW(&message as *const MSG);
            }
            true
        } else {
            false
        }
    }
}


fn main() {
    let mut window = create_window("my_window", "Simple Ray Tracer").unwrap();

    loop {
        if !handle_message(&mut window) {
            break;
        }
    }
}
