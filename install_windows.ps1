<#
.SYNOPSIS
    BACTolomew installer for Windows 10 (2004 or newer) and Windows 11, through WSL2.

.DESCRIPTION
    BACTolomew is a Linux program. On Windows it runs inside WSL2 (the Windows Subsystem
    for Linux), with its window shown by WSLg and its tools run by Docker. This script:

      1. checks the Windows version;
      2. checks WSL2, and installs it when missing (asks for administrator rights only for
         that step; Windows then needs a restart, after which you run this script again);
      3. checks or installs the Linux distribution (Ubuntu 24.04 by default);
      4. checks WSLg, which shows Linux windows on the Windows desktop;
      5. installs the few Linux packages the installer and the window need (git, Python's
         venv, the Qt system libraries), as root inside the distribution;
      6. checks Docker (not for -Viewer): Docker Desktop with its WSL integration, or
         Docker Engine installed inside the distribution;
      7. runs BACTolomew's own install.sh inside the distribution;
      8. creates a Start Menu shortcut (and a desktop shortcut with -DesktopShortcut).

    Every step checks first and does nothing when it is already done, so the script is
    safe to run again (after a restart, after installing Docker, or to repair).

.PARAMETER Viewer
    Install the window only, to open runs that someone else made: no Docker, tools image
    or databases (about 1 GB). `bactolomew setup` adds the pipeline later.

.PARAMETER Distro
    The WSL distribution to use (default: Ubuntu-24.04, or your existing default
    distribution after asking). It is installed when missing.

.PARAMETER Dir
    The install folder inside the distribution (default: ~/BACTolomew). Keep it in the
    Linux file system: folders under /mnt/c are much slower.

.PARAMETER Ref
    The BACTolomew release tag or branch to install (default: main).

.PARAMETER Docker
    How to get Docker: Ask (default), Desktop (Docker Desktop for Windows), Engine (Docker
    Engine inside the distribution) or Skip.

.PARAMETER DesktopShortcut
    Also put a BACTolomew shortcut on the desktop.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File .\install_windows.ps1

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File .\install_windows.ps1 -Viewer -DesktopShortcut
#>
[CmdletBinding()]
param(
    [switch]$Viewer,
    [string]$Distro = "",
    [string]$Dir = "",
    [string]$Ref = "main",
    [ValidateSet("Ask", "Desktop", "Engine", "Skip")]
    [string]$Docker = "Ask",
    [switch]$DesktopShortcut
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = "Stop"

$DefaultDistro = "Ubuntu-24.04"
$RepoRaw = "https://raw.githubusercontent.com/liviurotiul/BACTolomew"
$MinBuild = 19041           # Windows 10 version 2004
$MinWslgBuildWin10 = 19044  # WSLg on Windows 10 needs 21H2 or newer and the Store WSL
$SetupDiskGB = 150          # databases during setup (about 125 GB) plus the tools image (about 22 GB)
$RecommendedRamGB = 32

# --------------------------------------------------------------------------- #
# Messages
# --------------------------------------------------------------------------- #
function Say([string]$Text)  { Write-Host "[install] $Text" -ForegroundColor Cyan }
function Warn([string]$Text) { Write-Host "[install] WARNING: $Text" -ForegroundColor Yellow }
function Fail([string]$Text) {
    Write-Host "[install] ERROR: $Text" -ForegroundColor Red
    exit 1
}
function Ask-YesNo([string]$Question, [bool]$Default = $true) {
    $hint = "[y/N]"
    if ($Default) { $hint = "[Y/n]" }
    $answer = Read-Host "$Question $hint"
    if ([string]::IsNullOrWhiteSpace($answer)) { return $Default }
    return $answer.Trim().ToLower().StartsWith("y")
}

# --------------------------------------------------------------------------- #
# wsl.exe
# --------------------------------------------------------------------------- #
# A 32-bit PowerShell sees System32 redirected and finds no wsl.exe there.
$WslExe = Join-Path $env:SystemRoot "System32\wsl.exe"
$sysnative = Join-Path $env:SystemRoot "Sysnative\wsl.exe"
if (Test-Path $sysnative) { $WslExe = $sysnative }

# wsl.exe writes its own messages as UTF-16; WSL_UTF8 asks for UTF-8, and the NUL
# characters an older WSL still sends are dropped.
$env:WSL_UTF8 = "1"
try { [Console]::OutputEncoding = [Text.Encoding]::UTF8 } catch { }

# Runs wsl.exe with the given arguments; returns its output lines (NULs removed).
# The exit code is in $LASTEXITCODE afterwards.
function Wsl-Output([string[]]$Arguments) {
    if (-not (Test-Path $WslExe)) { $global:LASTEXITCODE = 1; return @() }
    # Windows PowerShell turns stderr lines into errors, which "Stop" would make fatal.
    $ErrorActionPreference = "Continue"
    $lines = & $WslExe @Arguments 2>&1
    $code = $LASTEXITCODE
    $clean = @($lines | ForEach-Object { ($_.ToString() -replace "`0", "").TrimEnd() } |
               Where-Object { $_ -ne "" })
    $global:LASTEXITCODE = $code
    return $clean
}

# A bash script for the distribution goes through base64 and a temporary file, so no
# quoting reaches wsl.exe and Windows line endings never reach bash.
function Distro-Command([string]$Script) {
    $unix = $Script -replace "`r`n", "`n"
    $b64 = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($unix))
    return "f=`$(mktemp); printf %s $b64 | base64 -d > `$f; bash `$f; rc=`$?; rm -f `$f; exit `$rc"
}

# Runs a bash script in the distribution and returns its output lines.
function Distro-Output([string]$Script, [switch]$AsRoot) {
    $a = @("-d", $script:Distro)
    if ($AsRoot) { $a += @("-u", "root") }
    $a += @("--exec", "bash", "-c", (Distro-Command $Script))
    return Wsl-Output $a
}

# Runs a bash script in the distribution in this console (its output shown live, its
# prompts answered here) and returns its exit code.
function Distro-Run([string]$Script, [switch]$AsRoot) {
    $argLine = "-d $script:Distro"
    if ($AsRoot) { $argLine += " -u root" }
    $argLine += " --exec bash -c `"" + (Distro-Command $Script) + "`""
    $p = Start-Process -FilePath $WslExe -ArgumentList $argLine -NoNewWindow -Wait -PassThru
    return $p.ExitCode
}

function Test-Admin {
    $id = [Security.Principal.WindowsIdentity]::GetCurrent()
    return (New-Object Security.Principal.WindowsPrincipal($id)).IsInRole(
        [Security.Principal.WindowsBuiltInRole]::Administrator)
}

# Runs one program elevated (a UAC prompt), waits, and returns its exit code.
function Run-Elevated([string]$File, [string]$Arguments) {
    if (Test-Admin) {
        $p = Start-Process -FilePath $File -ArgumentList $Arguments -NoNewWindow -Wait -PassThru
    } else {
        Say "Windows asks for administrator rights for: $File $Arguments"
        try {
            $p = Start-Process -FilePath $File -ArgumentList $Arguments -Verb RunAs -Wait -PassThru
        } catch {
            Fail "Administrator rights were not granted. Run this script again and accept the prompt, or run in an administrator PowerShell:  $File $Arguments"
        }
    }
    return $p.ExitCode
}

function Test-RebootPending {
    return (Test-Path "HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Component Based Servicing\RebootPending") -or
           (Test-Path "HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\WindowsUpdate\Auto Update\RebootRequired")
}

function Exit-ForRestart([string]$What) {
    Write-Host ""
    Say "$What"
    Say "Restart Windows now, then run this script again with the same options: it continues from here."
    exit 0
}

# --------------------------------------------------------------------------- #
# 1. Windows
# --------------------------------------------------------------------------- #
function Check-Windows {
    $build = [Environment]::OSVersion.Version.Build
    $name = "Windows 10"
    if ($build -ge 22000) { $name = "Windows 11" }
    $display = ""
    try {
        $cv = Get-ItemProperty "HKLM:\SOFTWARE\Microsoft\Windows NT\CurrentVersion"
        if ($cv.PSObject.Properties["DisplayVersion"]) { $display = $cv.DisplayVersion }
        elseif ($cv.PSObject.Properties["ReleaseId"]) { $display = $cv.ReleaseId }
    } catch { }
    if ($build -lt $MinBuild) {
        Fail "This is $name build $build. BACTolomew needs WSL2, which needs Windows 10 version 2004 (build $MinBuild) or newer, or Windows 11. Update Windows first (Settings > Windows Update)."
    }
    Say "Windows: $name $display (build $build)"
    $script:WinBuild = $build

    $arch = $env:PROCESSOR_ARCHITECTURE
    if ($env:PROCESSOR_ARCHITEW6432) { $arch = $env:PROCESSOR_ARCHITEW6432 }
    if ($arch -eq "ARM64") {
        Warn "This is a Windows on ARM PC. BACTolomew is tested on x86_64 only, and the published tools image is x86_64 only: on ARM, setup builds the tools image itself (several hours, without fastGEAR). A -Viewer install is the safer choice here."
    } elseif ($arch -ne "AMD64") {
        Fail "This is a $arch Windows. BACTolomew needs 64-bit Windows (x64)."
    }
}

# --------------------------------------------------------------------------- #
# 2. WSL2
# --------------------------------------------------------------------------- #
function Ensure-Wsl {
    $null = Wsl-Output @("--status")
    if ($LASTEXITCODE -ne 0) {
        if (Test-RebootPending) {
            Exit-ForRestart "Windows has a restart pending (WSL may have been installed already)."
        }
        Say "WSL is not installed (or not enabled). Installing it; this needs administrator rights once."
        $rc = Run-Elevated $WslExe "--install --no-distribution"
        if ($rc -ne 0) {
            Warn "wsl --install --no-distribution returned $rc. If this Windows is too old for that option, run in an administrator PowerShell:"
            Write-Host "    dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart"
            Write-Host "    dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart"
            Write-Host "  If it says virtualization is off, turn on Intel VT-x or AMD-V (SVM) in the PC's BIOS/UEFI settings."
        }
        Exit-ForRestart "WSL was installed (Windows features were enabled)."
    }

    # The Store version of WSL (it has `wsl --version`) brings WSLg and systemd support.
    $ver = @(Wsl-Output @("--version"))
    if ($LASTEXITCODE -ne 0) {
        Say "This is the older built-in WSL. Updating to the current WSL (it adds WSLg, which shows Linux windows)."
        $rc = Run-Elevated $WslExe "--update"
        $ver = @(Wsl-Output @("--version"))
        if ($LASTEXITCODE -ne 0) {
            Exit-ForRestart "wsl --update finished (exit code $rc), but the new WSL is not active yet."
        }
    }
    $wslVersion = ""
    if ($ver.Count -gt 0 -and $ver[0] -match "(\d+\.\d+\.\d+)") { $wslVersion = $Matches[1] }
    Say "WSL: version $wslVersion"
    if ($wslVersion -and ([version]$wslVersion -lt [version]"2.0.0")) {
        Warn "WSL $wslVersion is old. Update it with:  wsl --update"
    }
    # New distributions are WSL 2 (WSL 1 cannot run Docker).
    $null = Wsl-Output @("--set-default-version", "2")
}

# --------------------------------------------------------------------------- #
# 3. The distribution
# --------------------------------------------------------------------------- #
# Installed distributions: @{ Name; State; Version; Default }.
function Get-Distros {
    $rows = @()
    foreach ($line in (Wsl-Output @("--list", "--verbose"))) {
        $isDefault = $line.TrimStart().StartsWith("*")
        $parts = @(($line.TrimStart(" ", "*").Trim()) -split "\s+")
        if ($parts.Count -lt 3 -or $parts[-1] -notmatch "^\d$") { continue }   # the header line
        $rows += [pscustomobject]@{
            Name = $parts[0]; State = $parts[1]; Version = [int]$parts[-1]; Default = $isDefault
        }
    }
    return $rows
}

function Ensure-Distro {
    $distros = @(Get-Distros | Where-Object { $_.Name -notlike "docker-desktop*" })
    if (-not $script:Distro) {
        $preferred = $distros | Where-Object { $_.Name -eq $DefaultDistro }
        $default = $distros | Where-Object { $_.Default }
        if ($preferred) {
            $script:Distro = $DefaultDistro
        } elseif ($default) {
            $d = @($default)[0].Name
            Say "Your default WSL distribution is $d."
            Say "BACTolomew needs a glibc 2.34 or newer distribution for its window (Ubuntu 22.04 or 24.04, Debian 12, Fedora, openSUSE Tumbleweed or Leap 15.6)."
            if (Ask-YesNo "Install BACTolomew into $d? (No installs $DefaultDistro beside it)" $true) {
                $script:Distro = $d
            } else {
                $script:Distro = $DefaultDistro
            }
        } else {
            $script:Distro = $DefaultDistro
        }
    }

    $found = $distros | Where-Object { $_.Name -eq $script:Distro }
    if (-not $found) {
        $online = @(Wsl-Output @("--list", "--online"))
        $known = $online | Where-Object { ($_ -split "\s+")[0] -eq $script:Distro }
        if (-not $known) {
            Write-Host ($online -join "`n")
            Fail "'$($script:Distro)' is neither installed nor in the list above. Pick a name with -Distro."
        }
        Say "Installing the Linux distribution $($script:Distro) (a download of about 0.5 GB)."
        Say "It then asks for a new Linux user name and password: remember the password, Linux asks for it for administrator tasks."
        Say "When you see a Linux prompt (ending in `$), type  exit  and press Enter to come back here."
        $p = Start-Process -FilePath $WslExe -ArgumentList "--install -d $($script:Distro)" -NoNewWindow -Wait -PassThru
        $found = Get-Distros | Where-Object { $_.Name -eq $script:Distro }
        if (-not $found) {
            if (Test-RebootPending) { Exit-ForRestart "The distribution needs a Windows restart to finish installing." }
            Fail "$($script:Distro) did not install (wsl --install exit code $($p.ExitCode)). Run  wsl --install -d $($script:Distro)  by hand to see why, then run this script again."
        }
    }
    $found = @($found)[0]
    if ($found.Version -ne 2) {
        Say "$($script:Distro) runs on WSL 1; converting it to WSL 2 (a few minutes)."
        $null = Wsl-Output @("--set-version", $script:Distro, "2")
        if ($LASTEXITCODE -ne 0) { Fail "Could not convert $($script:Distro) to WSL 2. Run:  wsl --set-version $($script:Distro) 2" }
    }

    $info = @(Distro-Output 'printf "%s|%s|%s\n" "$(id -un)" "$(id -u)" "$HOME"')
    if ($LASTEXITCODE -ne 0 -or $info.Count -eq 0) {
        Fail "Could not start $($script:Distro). Try  wsl -d $($script:Distro)  by hand; if it never finished its first start, finish it there, then run this script again."
    }
    $fields = $info[-1] -split "\|"
    $script:LinuxUser = $fields[0]
    $script:LinuxHome = $fields[2]
    Say "Distribution: $($script:Distro) (WSL 2), Linux user $($script:LinuxUser)"
    if ($fields[1] -eq "0") {
        Warn "The distribution's default user is root. BACTolomew works as root, but a normal user is safer: create one inside the distribution (adduser), then make it the default (the [user] section of /etc/wsl.conf)."
    }
}

# --------------------------------------------------------------------------- #
# 4. WSLg
# --------------------------------------------------------------------------- #
function Check-Wslg {
    $out = @(Distro-Output 'if [ -d /mnt/wslg ] && [ -n "${WAYLAND_DISPLAY:-}${DISPLAY:-}" ]; then echo yes; else echo no; fi')
    if ($out -contains "yes") {
        Say "WSLg: Linux windows can open on the Windows desktop"
        return
    }
    $need = "Update WSL in a PowerShell:  wsl --update  then  wsl --shutdown  and run this script again."
    if ($script:WinBuild -lt 22000 -and $script:WinBuild -lt $MinWslgBuildWin10) {
        $need = "WSLg needs Windows 10 21H2 (build $MinWslgBuildWin10) or Windows 11; update Windows, then run  wsl --update."
    }
    Warn "WSLg is not available, so the BACTolomew window cannot open. $need"
    Warn "(Also check that %USERPROFILE%\.wslconfig has no  guiApplications=false  line.)"
}

# --------------------------------------------------------------------------- #
# 5. Linux packages (as root inside the distribution; no password needed)
# --------------------------------------------------------------------------- #
$PackagesScript = @'
set -e
. /etc/os-release
like=" ${ID:-} ${ID_LIKE:-} "
case "$like" in
  *" debian "*|*" ubuntu "*)
    t64=0
    case "$ID" in
      ubuntu) dpkg --compare-versions "${VERSION_ID:-0}" ge 24.04 && t64=1 ;;
      debian) dpkg --compare-versions "${VERSION_ID:-0}" ge 13 && t64=1 ;;
    esac
    asound=libasound2; glib=libglib2.0-0
    [ "$t64" = 1 ] && { asound=libasound2t64; glib=libglib2.0-0t64; }
    pkgs="git curl ca-certificates python3 python3-venv"
    pkgs="$pkgs libnss3 libnspr4 libxkbcommon0 libxkbcommon-x11-0 libxcomposite1 libxdamage1
      libxrandr2 libxfixes3 libxtst6 libxkbfile1 libgbm1 $asound libxshmfence1 fontconfig libgl1 libegl1 $glib
      libdbus-1-3 libgssapi-krb5-2 libxcb-cursor0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1
      libxcb-render-util0 libxcb-render0 libxcb-shape0 libxcb-util1 libxcb-xkb1 fonts-dejavu-core"
    missing=""
    for p in $pkgs; do dpkg -s "$p" >/dev/null 2>&1 || missing="$missing $p"; done
    if [ -z "$missing" ]; then echo "Linux packages: all present"; exit 0; fi
    echo "Installing:$missing"
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -q
    apt-get install -y -q $missing ;;
  *" fedora "*|*" rhel "*|*" centos "*|*" rocky "*|*" almalinux "*)
    pkgs="git curl python3"
    [ "$ID" = fedora ] || pkgs="$pkgs python3.12"
    pkgs="$pkgs nss nspr libxkbcommon libxkbcommon-x11 libXcomposite libXdamage libXrandr
      libXfixes libXtst libxkbfile mesa-libgbm alsa-lib libxshmfence fontconfig mesa-libGL mesa-libEGL glib2
      dbus-libs krb5-libs xcb-util-cursor xcb-util-wm xcb-util-image xcb-util-keysyms xcb-util-renderutil
      xcb-util dejavu-sans-fonts"
    if rpm -q --whatprovides $pkgs >/dev/null 2>&1; then echo "Linux packages: all present"; exit 0; fi
    dnf install -y $pkgs ;;
  *suse*)
    py=python3
    case "$ID" in opensuse-leap|sles) py=python312 ;; esac
    pkgs="git curl $py"
    pkgs="$pkgs mozilla-nss mozilla-nspr libxkbcommon0 libxkbcommon-x11-0 libXcomposite1
      libXdamage1 libXrandr2 libXfixes3 libXtst6 libxkbfile1 libgbm1 libasound2 libxshmfence1 fontconfig
      libglvnd libglib-2_0-0 libgthread-2_0-0 libdbus-1-3 krb5 libxcb-cursor0 libxcb-icccm4 libxcb-image0
      libxcb-keysyms1 libxcb-render-util0 libxcb-render0 libxcb-shape0 libxcb-util1 libxcb-xkb1
      dejavu-fonts"
    if rpm -q --whatprovides $pkgs >/dev/null 2>&1; then echo "Linux packages: all present"; exit 0; fi
    zypper --non-interactive install $pkgs ;;
  *)
    echo "Unknown distribution ($ID): install git, curl and Python 3.10 to 3.14 with venv yourself;"
    echo "install.sh names any missing library." ;;
esac
'@

function Install-Packages {
    Say "Checking the Linux packages the installer and the window need (as root inside $($script:Distro))"
    $rc = Distro-Run $PackagesScript -AsRoot
    if ($rc -ne 0) {
        Fail "Installing the Linux packages failed (exit code $rc). Check the network inside WSL (for example  wsl -d $($script:Distro) -- curl -I https://github.com), then run this script again."
    }
}

# --------------------------------------------------------------------------- #
# 6. Docker
# --------------------------------------------------------------------------- #
function Test-DistroDocker {
    $out = @(Distro-Output 'e=$(mktemp); if docker info >/dev/null 2>$e; then echo OK; else echo FAIL; tail -n 1 $e; fi; rm -f $e; command -v docker >/dev/null && echo HAVE_CLI || echo NO_CLI')
    return ,$out
}

# Docker Desktop: installed, running, and its WSL settings (read only).
function Get-DockerDesktop {
    $exe = Join-Path $env:ProgramFiles "Docker\Docker\Docker Desktop.exe"
    $state = [pscustomobject]@{
        Installed = (Test-Path $exe); Exe = $exe
        Running = [bool](Get-Process -Name "Docker Desktop" -ErrorAction SilentlyContinue)
        WslEngine = $null; Integrated = $false
    }
    foreach ($f in @("settings-store.json", "settings.json")) {
        $path = Join-Path $env:APPDATA "Docker\$f"
        if (-not (Test-Path $path)) { continue }
        try { $s = Get-Content $path -Raw | ConvertFrom-Json } catch { continue }
        foreach ($p in $s.PSObject.Properties) {
            if ($p.Name -eq "wslEngineEnabled") { $state.WslEngine = [bool]$p.Value }
            if ($p.Name -eq "integratedWslDistros" -and @($p.Value) -contains $script:Distro) { $state.Integrated = $true }
            if ($p.Name -eq "enableIntegrationWithDefaultWslDistro" -and $p.Value) {
                $def = Get-Distros | Where-Object { $_.Default -and $_.Name -eq $script:Distro }
                if ($def) { $state.Integrated = $true }
            }
        }
        break
    }
    return $state
}

$DockerEngineScript = @'
set -e
. /etc/os-release
like=" ${ID:-} ${ID_LIKE:-} "
# Docker Engine needs systemd in WSL (WSL 0.67.6 or newer reads [boot] systemd=true).
if [ "$(cat /proc/1/comm 2>/dev/null)" != systemd ]; then
  touch /etc/wsl.conf
  if grep -qE '^[[:space:]]*systemd[[:space:]]*=' /etc/wsl.conf; then
    sed -i -E 's/^[[:space:]]*systemd[[:space:]]*=.*/systemd=true/' /etc/wsl.conf
  elif grep -q '^\[boot\]' /etc/wsl.conf; then
    sed -i '/^\[boot\]/a systemd=true' /etc/wsl.conf
  else
    printf '\n[boot]\nsystemd=true\n' >> /etc/wsl.conf
  fi
  echo "systemd turned on in /etc/wsl.conf; the distribution restarts."
  exit 10
fi
if ! command -v dockerd >/dev/null 2>&1; then
  case "$like" in
    *suse*) zypper --non-interactive install docker ;;
    *)
      tmp="$(mktemp -d)"
      curl -fsSL -o "$tmp/get-docker.sh" https://get.docker.com
      sh "$tmp/get-docker.sh"
      rm -rf "$tmp" ;;
  esac
fi
systemctl enable --now docker
getent group docker >/dev/null || groupadd docker
if [ -n "$TARGET_USER" ] && [ "$TARGET_USER" != root ]; then usermod -aG docker "$TARGET_USER"; fi
echo "Docker Engine installed and started."
'@

function Install-DockerEngine {
    $dd = Get-DockerDesktop
    if ($dd.Installed) {
        Warn "Docker Desktop is installed too. Do not turn on its WSL integration for $($script:Distro) while Docker Engine runs inside it: the two conflict."
    }
    $userLine = "TARGET_USER='" + ($script:LinuxUser -replace "'", "") + "'`n"
    for ($i = 0; $i -lt 2; $i++) {
        $rc = Distro-Run ($userLine + $DockerEngineScript) -AsRoot
        if ($rc -ne 10) { break }
        $null = Wsl-Output @("--terminate", $script:Distro)
    }
    if ($rc -eq 10) {
        Fail "systemd did not start in $($script:Distro) after a restart; Docker Engine needs it. Update WSL (wsl --update, then wsl --shutdown) and run this script again, or use Docker Desktop (-Docker Desktop)."
    }
    if ($rc -ne 0) {
        Fail "Installing Docker Engine inside $($script:Distro) failed (exit code $rc). See the messages above; https://docs.docker.com/engine/install/ has the steps per distribution."
    }
    # A new start of the distribution picks up the docker group.
    Say "Restarting $($script:Distro) so that $($script:LinuxUser) can use Docker (this closes its open terminals)."
    $null = Wsl-Output @("--terminate", $script:Distro)
}

function Install-DockerDesktop {
    $winget = Get-Command winget.exe -ErrorAction SilentlyContinue
    if (-not $winget) {
        Warn "winget is not available. Download Docker Desktop from https://docs.docker.com/desktop/setup/install/windows-install/ and install it."
        return
    }
    Say "Installing Docker Desktop with winget (it may ask for administrator rights)."
    & $winget.Source install -e --id Docker.DockerDesktop --accept-package-agreements --accept-source-agreements
    if ($LASTEXITCODE -ne 0) {
        Warn "winget returned $LASTEXITCODE. Install Docker Desktop from https://docs.docker.com/desktop/setup/install/windows-install/ instead."
    }
}

function Show-DockerDesktopSteps {
    Write-Host "  In Docker Desktop:"
    Write-Host "    1. start it once and finish its first-run steps (Docker's terms; it is free for personal use,"
    Write-Host "       education and small companies, larger companies need a paid subscription);"
    Write-Host "    2. Settings > General: 'Use the WSL 2 based engine' on;"
    Write-Host "    3. Settings > Resources > WSL integration: turn on $($script:Distro), then Apply & restart."
    Write-Host "  Then run this script again (or  ./bactolomew setup  inside $($script:Distro))."
}

function Ensure-Docker {
    if ($Viewer) {
        Say "Viewer install: Docker is not needed (running predictions later needs it; run this script again without -Viewer)."
        return
    }
    if ($Docker -eq "Skip") {
        Say "Skipping the Docker check (-Docker Skip)."
        return
    }
    $out = @(Test-DistroDocker)
    if ($out -contains "OK") {
        Say "Docker: usable inside $($script:Distro)"
        return
    }
    $dd = Get-DockerDesktop
    $hasCli = $out -contains "HAVE_CLI"
    $err = ($out | Where-Object { $_ -notin @("OK", "FAIL", "HAVE_CLI", "NO_CLI") }) -join " "

    if ($dd.Installed -and $Docker -ne "Engine") {
        Warn "Docker Desktop is installed but Docker does not work inside $($script:Distro)."
        if (-not $dd.Running) { Write-Host "  Docker Desktop is not running: start it from the Start Menu and wait until it says it is running." }
        if ($dd.WslEngine -eq $false) { Write-Host "  'Use the WSL 2 based engine' is off." }
        if (-not $dd.Integrated) { Write-Host "  The WSL integration for $($script:Distro) looks off." }
        if ($err) { Write-Host "  docker info said: $err" }
        Show-DockerDesktopSteps
        return
    }
    if ($hasCli -and $Docker -eq "Ask") {
        Warn "Docker is installed inside $($script:Distro) but not usable: $err"
        Write-Host "  If its service is stopped:  sudo systemctl start docker"
        Write-Host "  If permission was denied:  sudo usermod -aG docker `$USER  then  wsl --terminate $($script:Distro)"
        return
    }

    $choice = $Docker
    if ($choice -eq "Ask") {
        Write-Host ""
        Say "BACTolomew runs its tools with Docker. Two ways to get it on Windows:"
        Write-Host "  1. Docker Desktop for Windows (recommended: the simplest, with a settings window and updates;"
        Write-Host "     free for personal use, education and small companies)"
        Write-Host "  2. Docker Engine inside $($script:Distro) (free for everyone, command line only; turns on systemd in the distribution)"
        Write-Host "  3. Not now (BACTolomew installs; analyses need Docker later)"
        $a = Read-Host "Choose 1, 2 or 3 [1]"
        switch ($a.Trim()) {
            "2" { $choice = "Engine" }
            "3" { $choice = "Skip" }
            default { $choice = "Desktop" }
        }
    }
    switch ($choice) {
        "Desktop" {
            if (-not $dd.Installed) { Install-DockerDesktop }
            Show-DockerDesktopSteps
        }
        "Engine" {
            Install-DockerEngine
            $out = @(Test-DistroDocker)
            if ($out -contains "OK") { Say "Docker: usable inside $($script:Distro)" }
            else { Warn "Docker Engine was installed but does not answer yet. Run this script again; if that fails, run  wsl --shutdown  first." }
        }
        default { Say "Docker skipped: install it later and run this script again." }
    }
}

# --------------------------------------------------------------------------- #
# 7. install.sh inside the distribution
# --------------------------------------------------------------------------- #
function Quote-Bash([string]$s) { return "'" + ($s -replace "'", "'\''") + "'" }

$InstallScript = @'
set -e
case "$DIR" in "~") DIR="$HOME" ;; "~/"*) DIR="$HOME/${DIR#\~/}" ;; esac
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
url="$RAW/$REF/install.sh"
echo "Downloading $url"
curl -fsSL -o "$tmp/install.sh" "$url"
set -- --ref "$REF"
[ -n "$DIR" ] && set -- "$@" --dir "$DIR"
[ "$VIEWER" = 1 ] && set -- "$@" --viewer
cd "$HOME"
bash "$tmp/install.sh" "$@"
'@

function Run-Installer {
    if ($Dir -like "/mnt/*") {
        Warn "$Dir is on a Windows drive: analyses there run many times slower than in the Linux file system (for example ~/BACTolomew)."
        if (-not (Ask-YesNo "Install there anyway?" $false)) { Fail "Run the script again with another -Dir (or none)." }
    }
    $vars = "RAW=" + (Quote-Bash $RepoRaw) + "`nREF=" + (Quote-Bash $Ref) + "`nDIR=" + (Quote-Bash $Dir) +
            "`nVIEWER=" + [int][bool]$Viewer + "`n"
    Say "Running BACTolomew's install.sh ($Ref) inside $($script:Distro)"
    $rc = Distro-Run ($vars + $InstallScript)
    if ($rc -ne 0) { Fail "install.sh failed (exit code $rc); its messages are above. Fix what it names, then run this script again." }

    $find = "DIR=" + (Quote-Bash $Dir) + "`n" + @'
case "$DIR" in "") DIR="$HOME/BACTolomew" ;; "~") DIR="$HOME" ;; "~/"*) DIR="$HOME/${DIR#\~/}" ;; esac
[ -x "$DIR/bactolomew" ] && (cd "$DIR" && echo "$PWD/bactolomew")
'@
    $launcher = @(Distro-Output $find)
    if ($launcher.Count -eq 0) { Fail "install.sh finished, but no bactolomew launcher was found in the install folder." }
    $script:Launcher = $launcher[-1]
}

# --------------------------------------------------------------------------- #
# 8. Shortcuts
# --------------------------------------------------------------------------- #
function New-BactolomewShortcut([string]$Folder) {
    $path = Join-Path $Folder "BACTolomew.lnk"
    $shell = New-Object -ComObject WScript.Shell
    $lnk = $shell.CreateShortcut($path)
    $lnk.TargetPath = Join-Path $env:SystemRoot "System32\wsl.exe"
    $lnk.Arguments = "-d $($script:Distro) --cd ~ --exec `"$($script:Launcher)`""
    $lnk.WorkingDirectory = $env:USERPROFILE
    $lnk.Description = "BACTolomew (runs in WSL: $($script:Distro))"
    $lnk.WindowStyle = 7   # the console window starts minimized
    $lnk.Save()
    Say "Shortcut: $path"
}

function Create-Shortcuts {
    New-BactolomewShortcut ([Environment]::GetFolderPath("Programs"))
    if ($DesktopShortcut) { New-BactolomewShortcut ([Environment]::GetFolderPath("Desktop")) }
}

# --------------------------------------------------------------------------- #
# 9. Memory and disk
# --------------------------------------------------------------------------- #
function Check-Resources {
    $r = @(Distro-Output @'
read -r _ kb _ < /proc/meminfo; echo $(( kb / 1048576 ))
set -- $(df -Pk "$HOME" | tail -n 1); echo $(( $4 / 1048576 ))
'@)
    $hostRam = [int]((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB)
    if ($r.Count -ge 2) {
        $wslRam = [int]$r[0]
        Say "Memory: WSL sees $wslRam GB of the PC's $hostRam GB"
        if (-not $Viewer -and $wslRam -lt ($RecommendedRamGB - 2)) {
            if ($hostRam -gt ($wslRam + 4)) {
                $give = [Math]::Min($RecommendedRamGB, $hostRam - 4)
                Warn "BACTolomew analyses want about $RecommendedRamGB GB. WSL gets half of the PC's memory by default; to give it more, put these lines in $env:USERPROFILE\.wslconfig and run  wsl --shutdown :"
                Write-Host "    [wsl2]"
                Write-Host "    memory=${give}GB"
            } else {
                Warn "BACTolomew analyses want about $RecommendedRamGB GB; with $hostRam GB, CheckM2 runs in its low-memory mode and big runs are slower."
            }
        }
    }

    if ($Viewer) { return }
    # The distribution's virtual disk (and Docker Desktop's) grow on a Windows drive: check that drive.
    $base = $null
    try {
        foreach ($key in (Get-ChildItem "HKCU:\Software\Microsoft\Windows\CurrentVersion\Lxss")) {
            $p = Get-ItemProperty $key.PSPath
            if ($p.PSObject.Properties["DistributionName"] -and $p.DistributionName -eq $script:Distro) { $base = $p.BasePath }
        }
    } catch {
        $base = $null   # no registry entry: the disk check is skipped
    }
    if ($base) {
        $root = [IO.Path]::GetPathRoot(($base -replace "^\\\\\?\\", ""))
        $drive = Get-PSDrive -Name $root.Substring(0, 1) -ErrorAction SilentlyContinue
        if ($drive) {
            $free = [int]($drive.Free / 1GB)
            Say "Disk: $free GB free on $root (where the $($script:Distro) virtual disk grows)"
            if ($free -lt $SetupDiskGB) {
                Warn "First-time setup needs about $SetupDiskGB GB there (databases about 125 GB while unpacking, the tools image about 22 GB). Free space on $root, or move the distribution:  wsl --manage $($script:Distro) --move <folder on another drive>"
            }
        }
    }
}

# --------------------------------------------------------------------------- #
function Main {
    Write-Host ""
    Say "BACTolomew for Windows (through WSL2)"
    Check-Windows
    Ensure-Wsl
    Ensure-Distro
    Check-Wslg
    Install-Packages
    Ensure-Docker
    Run-Installer
    Create-Shortcuts
    Check-Resources

    Write-Host ""
    Say "Done. Start BACTolomew from the Start Menu (BACTolomew), or in a PowerShell:"
    Write-Host "    wsl -d $($script:Distro) --cd ~ --exec $($script:Launcher)"
    Say "Keep genomes and runs inside Linux (for example \\wsl.localhost\$($script:Distro)$($script:LinuxHome) in File Explorer): files under C:\ (/mnt/c) work but are much slower."
    if ($Viewer) {
        Say "Viewer install: open shared runs with Import Run... in the run list."
    } else {
        Say "The first start gets the tools image and the databases (about 150 GB, several hours)."
    }
}

Main
