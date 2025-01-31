Graphics = love.graphics
Audio = love.audio
Keyboard = love.keyboard
Window = love.window
Mouse = love.mouse

sqrt = math.sqrt
abs = math.abs
cos = math.cos
sin = math.sin
exp = math.exp
ceil = math.ceil
floor = math.floor
min = math.min
max = math.max
deg_to_rad = math.pi / 180
rad_to_deg = 180 / math.pi

Clamp = function (x, u, v)
	return max(u, min(x, v))
end

Linear = function (x, x0, y0, x1, y1)
	local f = (x - x0) / (x1 - x0)
	return f * y1 + (1 - f) * y0
end

function MakeClass (m, base)
	m.__index = base or m
	setmetatable(m, {__index = base})
	return m
end

function Create (what)
	local m = setmetatable({}, {__index = what})
	if m.OnBuild then
		m:OnBuild()
	end
	return m
end

Grid = MakeClass {
	Copy = function (m, other)
		m.x_size = other.x_size
		m.y_size = other.y_size
		
		for x = 1, m.x_size do
			for y = 1, m.y_size do
				m [x] [y] = other [x] [y]
			end
		end
	end,
	
	Fill = function (m, v)
		for x = 1, m.x_size do
			for y = 1, m.y_size do
				m [x] [y] = v
			end
		end
	end,
	
	Bilerp = function (m, x, y)
		x = Clamp(x, 1, m.x_size - 1)
		y = Clamp(y, 1, m.y_size - 1)
		
		local xi = floor(x)
		local yi = floor(y)
		local xf = x - xi
		local yf = y - yi
		local uf = 1 - xf
		local vf = 1 - yf
		
		return
			xf * yf * m [xi + 1] [yi + 1] +
			uf * yf * m [xi + 0] [yi + 1] +
			xf * vf * m [xi + 1] [yi + 0] +
			uf * vf * m [xi + 0] [yi + 0]
	end,
	
	Resize = function (m, x_size, y_size)
		x_size = x_size or 1
		y_size = y_size or 1
		
		for x = 1, x_size do
			m [x] = {}
			for y = 1, y_size do
				m [x] [y] = 0
			end
		end
		
		m.x_size = x_size
		m.y_size = y_size
		return m
	end,
	
	New = function (x, y)
		local m = setmetatable({}, Grid)
		m:Resize(x, y)
		return m
	end,
}

function DifferentSigns (u, v)
	return (u <= 0 and 0 <= v) or (v <= 0 and 0 <= u)
end

function love.load ()
	love.app = Create(App)
end

function love.update (dt)
	love.app:OnTick()
end

function love.draw ()
	love.app:OnDraw()
end

State = MakeClass {
	OnBuild = function (m)
		m.timer = 1
	end,
	
	OnFirstTick = function (m) end,
	OnTick = function (m) end,
	OnDraw = function (m) Graphics.print("Hello default state!") end,
}

App = MakeClass {
	OnBuild = function (m)
		m.data = {}
		local d = m.data
		d.version = "ver. xyz"
		m.state_lib = MakeStateLib(d)
		
		-- m:GoToState(State)
		m:GoToState(m.state_lib.Simulate)
		-- m:GoToState(m.state_lib.SatTest)
		-- print(m.state_lib.Playing.OnDraw, State.OnDraw)
	end,
	
	GoToState = function (m, what_state)
		m.pending_state = m.state_lib [what_state] or what_state
	end,
	
	OnTick = function (m)
		if m.pending_state then
			m.current_state = Create(m.pending_state)
			m.pending_state = nil
		end
		
		local state = m.current_state

		if not state then
			return
		end

		if 1 == state.timer then
			state:OnFirstTick()
		end

		state:OnTick()
		state.timer = 1 + state.timer;
	end,
	
	OnDraw = function (m)
		local state = m.current_state
		
		if not state then
			return
		end
		
		state:OnDraw()
	end
}

function MakeSimulationLib (d)
	local lib = {}
	
	function lib.Divergence (hor, ver, x, y)
		return hor [x + 1] [y] - hor [x] [y] + ver [x] [y + 1] - ver [x] [y]
	end
	
	function lib.StaggeredAverage (grid, x, y)
		return 0.25 * (grid [x] [y] + grid [1 + x] [y] + grid [x] [1 + y] + grid [x + 1] [y + 1])
	end
	
	function lib.Div0ConditionXyp (hor, ver, pressure, free, h, density, dt, over_relax, iterations)
		local x_size = min(hor.x_size, ver.x_size, pressure.x_size, free.x_size)
		local y_size = min(hor.y_size, ver.y_size, pressure.y_size, free.y_size)
		
		-- Gauß-Seidel iteration
		for k = 1, iterations do
			
			-- Calculate divergence by counting obstacles. Ignore bordering cells as though they
			-- were obstacles. They are updated by ExtrapolateXy instead.
			for x = 2, x_size - 1 do
				for y = 2, y_size - 1 do
					-- Count how many neighbouring blocks are free.
					local free_l = free [x - 1] [y]
					local free_r = free [x + 1] [y]
					local free_t = free [x] [y - 1]
					local free_b = free [x] [y + 1]
					local free_count = free_l + free_r + free_t + free_b
					
					-- If the cell is trapped, there is no flow. Otherwise continue.
					if 0 == free_count then
						goto continue
					end
					
					local div = lib.Divergence(hor, ver, x, y)
					local pressure_base = over_relax * div / free_count
					local pressure_real = density * h * pressure_base / dt
					
					-- Update horizontal and vertical staggered velocities such that divergence will
					-- converge to 0 in all cells. This is a linear equation system being solved
					-- during the loop iteratively.
					hor [x] [y] = hor [x] [y] + free_l * pressure_base
					hor [x + 1] [y] = hor [x + 1] [y] - free_r * pressure_base
					ver [x] [y] = ver [x] [y] + free_t * pressure_base
					ver [x] [y + 1] = ver [x] [y + 1] - free_b * pressure_base
					
					pressure [x] [y] = pressure_real
					
					::continue::
				end
			end
		end
	end
	
	function lib.ExtrapolateXy (hor, ver)
		for y = 1, hor.y_size do
			hor [1] [y] = hor [2] [y]
			hor [hor.x_size] [y] = hor [hor.x_size - 1] [y]
		end
		
		for x = 1, ver.x_size do
			ver [x] [1] = ver [x] [2]
			ver [x] [ver.y_size] = ver [x] [ver.y_size - 1]
		end
	end
	
	function lib.Gravity (grid, gravity)
		for x = 1, grid.x_size do
			for y = 1, grid.y_size do
				grid [x] [y] = gravity + grid [x] [y]
			end
		end
	end
	
	-- Advection of velocities.
	function lib.AdvectXy (hor, ver, new_hor, new_ver, h, dt)
		local x_size = min(hor.x_size, ver.x_size)
		local y_size = min(hor.y_size, ver.y_size)
		
		for x = 2, x_size - 1 do
			for y = 2, y_size - 1 do
				local x_vel = lib.StaggeredAverage(hor, x, y)
				local y_vel = lib.StaggeredAverage(ver, x, y)
				local xf = x + 0.0 - dt * x_vel
				local yf = y + 0.0 - dt * y_vel
				
				new_hor [x] [y] = hor:Bilerp(xf, yf)
				new_ver [x] [y] = ver:Bilerp(xf, yf)
			end
		end
		
		hor:Copy(new_hor)
		ver:Copy(new_ver)
	end
	
	-- Using Runge-Kutta method.
	function lib.AdvectXyRk (hor, ver, new_hor, new_ver, h, dt)
		local x_size = min(hor.x_size, ver.x_size)
		local y_size = min(hor.y_size, ver.y_size)
		
		for x = 2, x_size - 1 do
			for y = 2, y_size - 1 do
				local x_vel = lib.StaggeredAverage(hor, x, y)
				local y_vel = lib.StaggeredAverage(ver, x, y)
				
				-- This is the implicit Euler method. I'd like using Runge-Kutta.
				local xf = x - dt * x_vel
				local yf = y - dt * y_vel
				
				new_hor [x] [y] = hor:Bilerp(xf, yf)
				new_ver [x] [y] = ver:Bilerp(xf, yf)
			end
		end
		
		hor:Copy(new_hor)
		ver:Copy(new_ver)
	end
	
	-- Advection of a custom property by velocities.
	function lib.AdvectXyProp (hor, ver, prop, new_prop, h, dt)
		local x_size = min(hor.x_size, ver.x_size, prop.x_size)
		local y_size = min(hor.y_size, ver.y_size, prop.y_size)
		
		for x = 2, x_size - 1 do
			for y = 2, y_size - 1 do
				local x_vel = lib.StaggeredAverage(hor, x, y)
				local y_vel = lib.StaggeredAverage(ver, x, y)
				local xf = x + 0.0 - dt * x_vel
				local yf = y + 0.0 - dt * y_vel
				
				new_prop [x] [y] = prop:Bilerp(xf, yf)
				-- new_prop [x] [y] = prop [Clamp(floor(xf), 1, x_size)] [Clamp(floor(yf), 1, y_size)]
			end
		end
		
		prop:Copy(new_prop)
	end
	
	function lib.Draw (grid)
		for x = 1, grid.x_size do
			for y = 1, grid.y_size do
				local c = grid [x] [y]
				Graphics.setColor(1, 0, 0, Clamp(c, 0, 1))
				Graphics.points(x, y)
			end
		end
	end
	
	function lib.Timer ()
		return d.timer
	end
	
	-- Primitive
	function lib.MakeXySim (x_size, y_size)
		x_size = x_size or 256
		y_size = y_size or 128
		
		return MakeClass {
			OnBuild = function (m)
				
				-- Cursor info
				m.UpdateCursor = function ()
					m.old_cursor_x, m.old_cursor_y = m.cursor_x or 0, m.cursor_y or 0
					m.cursor_x, m.cursor_y = Mouse.getPosition()
					
					-- Dependent on outside things.
					local base_x = 32
					local base_y = 32
					local scale = 3
					
					m.cursor_x = (m.cursor_x - base_x) / scale
					m.cursor_y = (m.cursor_y - base_y) / scale
				end
				
				m.UpdateCursor()
			
				-- Create a bunch of grids to store the velocities and pressure in.
				m.x_size = x_size
				m.y_size = y_size
				
				local grids = { "x_vel", "y_vel", "new_x_vel", "new_y_vel", "pressure", "smoke", "new_smoke", "free" }
				
				for k, v in ipairs(grids) do
					print(v, m.x_size, m.y_size)
					m [v] = Grid.New(m.x_size, m.y_size)
				end
				
				m.free:Fill(1)
			end,
			
			OnTick = function (m)
				local h = 0.5
				local gravity = 0.005
				local density = 2
				local iterations = 5
				local over_relax = 1.8
				local dt = 4.0
				
				m.UpdateCursor()
				
				-- Create a pipe of in-flow.
				for k = 60, m.y_size - 60 do
					m.x_vel [1] [k] = 0.8
					m.x_vel [2] [k] = 0.8
					m.smoke [1] [k] = 0.9
					m.smoke [2] [k] = 0.9
				end
				
				local speed = 0.25
				local dx = m.cursor_x - m.old_cursor_x
				local dy = m.cursor_y - m.old_cursor_y
				local vx = speed * dx
				local vy = speed * dy
				
				m:CircleObstacle(floor(m.cursor_x), floor(m.cursor_y), 15, vx, vy)
				
				-- lib.Gravity(m.y_vel, gravity)
				m.pressure:Fill(0)
				lib.Div0ConditionXyp(m.x_vel, m.y_vel, m.pressure, m.free, h, density, dt, over_relax, iterations)
				lib.AdvectXyRk(m.x_vel, m.y_vel, m.new_x_vel, m.new_y_vel, h, dt)
				lib.AdvectXyProp(m.x_vel, m.y_vel, m.smoke, m.new_smoke, h, dt)
			end,
			
			CircleObstacle = function (m, cx, cy, radius, vx, vy)
				local x1 = floor(cx - radius)
				local y1 = floor(cy - radius)
				local x2 = floor(cx + radius)
				local y2 = floor(cy + radius)
				local radius_squared = radius * radius
				
				x1 = 1
				x2 = x_size
				y1 = 1
				y2 = y_size
				vx = vx or 0
				vy = vy or 0
				
				for x = x1, x2 do
					for y = y1, y2 do
						local dx = x - cx
						local dy = y - cy
						
						if dx * dx + dy * dy <= radius_squared then
							m.free [x] [y] = 0
							-- m.smoke [x] [y] = 1 + 0.5 * sin(lib.Timer() * 0.05)
							-- m.x_vel [x] [y] = m.x_vel [x] [y] + vx
							-- m.y_vel [x] [y] = m.y_vel [x] [y] + vy
							m.x_vel [x] [y] = vx
							m.y_vel [x] [y] = vy
						else
							m.free [x] [y] = 1
							-- m.smoke [x] [y] = 0
						end
					end
				end
			end,
			
			DrawComposite = function (m)
				function TempColor (f)
					f = Clamp(f, 0, 1)
					
					if f < 0.45 then
						return { 0, Linear(f, 0, 0, 0.45, 1), Linear(f, 0, 0.6, 0.45, 1) }
					elseif f < 0.5 then
						return { 0, 1, Linear(f, 0.45, 1, 0.5, 0) }
					elseif f < 0.6 then
						return { Linear(f, 0.5, 0, 0.6, 1), 1, 0 }
					end
					
					return { 1, Linear(f, 0.6, 1, 1, 0), 0 }
				end
				
				
				for x = 1, x_size do
					for y = 1, y_size do
						-- local f = Clamp((abs(m.x_vel [x] [y]) + abs(m.x_vel [x] [y])) * 0.5, 0, 1)
						-- local f = Clamp(abs(m.pressure [x] [y]) * 1, 0, 1)
						local f = Clamp(abs(m.smoke [x] [y]) * 1, 0, 1)
						local c = TempColor(f)
						
						if 0 < m.free [x] [y] then
							Graphics.setColor(c [1], c [2], c [3], 1)
						else
							Graphics.setColor(0.6, 0.3, 0.8030, 1)
						end
						
						Graphics.points(x, y)
					end
				end
			end
		}
	end
	
	return lib
end

function MakeStateLib (d)
	local lib = {}
	
	lib.Simulate = MakeClass({
		
		OnFirstTick = function (m)
			d.timer = 1
			m.x_size = 240
			m.y_size = 148
			m.simulation_lib = MakeSimulationLib(d)
			m.sim = Create(m.simulation_lib.MakeXySim(m.x_Size, m.y_size))
			m.image = Graphics.newCanvas(m.x_size, m.y_size)
		end,
		
		OnTick = function (m)
			if 200 < d.timer then
				m.sim:OnTick()
			end
			
			d.timer = 1 + d.timer
		end,
		
		OnDraw = function (m)
			local function Fill ()
				local canvas = Graphics.getCanvas()
				local x_size, y_size
				
				if canvas then
					x_size, y_size = canvas:getDimensions()
				else
					x_size, y_size = Window.getMode()
				end
				
				Graphics.rectangle("fill", 0, 0, x_size, y_size)
			end
			
			-- Copy contents of the simulation into the image.
			Graphics.setColor(0, 0, 0, 1)
			Graphics.setCanvas(m.image)
			Graphics.clear()
			Fill()
			m.sim:DrawComposite()
			Graphics.setCanvas()
			
			Graphics.clear()
			Graphics.setColor(0.07, 0.03, 0.3, 1)
			Fill()
			Graphics.setColor(1, 1, 1, 1)
			Graphics.setLineStyle("rough")
			Graphics.print("Simulation Running.")
			Graphics.print(m.timer or 0, 0, 16)
			Graphics.push()
			Graphics.translate(32, 90)
			Graphics.scale(3, 3)
			Graphics.draw(m.image, 0, 0)
			Graphics.setColor(0, 1, 1, 0.6)
			Graphics.rectangle("line", 0.5, 0.5, m.x_size, m.y_size)
			Graphics.setColor(1, 1, 1, 1)
			Graphics.pop()
		end
		
	}, State)
	
	return lib
end
