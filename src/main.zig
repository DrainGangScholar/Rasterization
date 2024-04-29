const std = @import("std");
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const WIDTH: usize = 640;
const HEIGHT: usize = 480;
const CHANNELS: usize = 3;
const stdout = std.io.getStdOut();
const Point = struct {
    x: i64,
    y: i64,
    pub fn x_u64(self: *const Point) u64 {
        return @intCast(self.x);
    }
    pub fn y_u64(self: *const Point) u64 {
        return @intCast(self.y);
    }
};
pub fn BresenhamLine(start: Point, end: Point, allocator: Allocator) !ArrayList(Point) {
    var points = ArrayList(Point).init(allocator);
    errdefer points.deinit();
    const dx: i64 = @intCast(@abs(start.x - end.x));
    const dy: i64 = @intCast(@abs(start.y - end.y));

    const start_greater: bool = start.x > end.x;
    var x = if (start_greater) end.x else start.x;
    var y = if (start_greater) end.y else start.y;
    const xend = if (start_greater) start.x else end.x;

    var d: i64 = 2 * dy - dx;
    const incr1: i64 = 2 * dy;
    const incr2: i64 = 2 * (dy - dx);

    try points.append(Point{
        .x = x,
        .y = y,
    });

    while (x < xend) {
        x += 1;
        if (d < 0) {
            d += incr1;
        } else {
            y += 1;
            d += incr2;
        }
        try points.append(Point{
            .x = x,
            .y = y,
        });
    }
    return points;
}
pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const start = Point{ .x = 0, .y = 0 };
    const end = Point{ .x = WIDTH, .y = HEIGHT };

    const points = try BresenhamLine(start, end, allocator);
    defer points.deinit();

    for (points.items) |point| {
        std.debug.print("{}\n", .{point});
    }
    const write = stdout.writer();
    const size = CHANNELS * WIDTH * HEIGHT;
    var buffer = try allocator.alloc(u8, size);
    defer allocator.free(buffer);
    try write.print("P3\n{} {}\n255\n", .{ WIDTH, HEIGHT });
    var index: usize = 0;
    for (0..HEIGHT) |i| {
        for (0..WIDTH) |j| {
            index = (i * WIDTH + j) * CHANNELS;
            buffer[index] = 255;
            buffer[index + 1] = 255;
            buffer[index + 2] = 255;
            //try write.print("{} {} {}\n", .{ buffer[index], buffer[index + 1], buffer[index + 2] });
        }
    }

    index = 0;
    for (points.items) |point| {
        index = (point.y_u64() * WIDTH + point.x_u64()) * CHANNELS;
        try write.print("Index {}", .{index});
        buffer[index] = 0;
        buffer[index + 1] = 0;
        buffer[index + 2] = 0;
    }

    for (0..HEIGHT) |i| {
        for (0..WIDTH) |j| {
            index = (i * WIDTH + j) * CHANNELS;
            try write.print("{} {} {}\n", .{ buffer[index], buffer[index + 1], buffer[index + 2] });
        }
    }
}
