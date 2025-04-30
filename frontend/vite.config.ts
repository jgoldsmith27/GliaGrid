import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react-swc'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  server: {
    open: false, // Prevent opening browser tab automatically
    port: 5173, // Explicitly define frontend port (optional but good practice)
    proxy: {
      // Proxy requests starting with /api to the backend
      '/api': {
        target: 'http://127.0.0.1:8000', // Explicitly use IPv4 loopback address
        changeOrigin: true,
        // No rewrite needed, backend expects /api
      }
    }
  }
})
